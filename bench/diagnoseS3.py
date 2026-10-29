#!/usr/bin/env python3
"""
Why /vsis3 says "No such file or directory" when the object is plainly there.

Splits the two possible causes, which look identical from GDAL:
  * the credentials or permissions are wrong        -> boto3 fails too
  * GDAL is configured wrong / needs ListBucket     -> boto3 succeeds, GDAL does not

ASF's in-region role typically grants GetObject but NOT ListBucket, and GDAL's stat can fall
back to a directory listing, which then 403s and is reported as a missing file.

    python diagnoseS3.py <granule-stem>
"""
import os
import sys


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    granule = sys.argv[1]

    import earthaccess
    earthaccess.login(strategy='netrc')
    c = earthaccess.get_s3_credentials(daac='ASF')
    env = {'AWS_ACCESS_KEY_ID': c['accessKeyId'],
           'AWS_SECRET_ACCESS_KEY': c['secretAccessKey'],
           'AWS_SESSION_TOKEN': c['sessionToken'],
           'AWS_REGION': 'us-west-2', 'AWS_DEFAULT_REGION': 'us-west-2'}
    os.environ.update(env)
    print('credentials: id %d, key %d, token %d chars'
          % (len(c['accessKeyId']), len(c['secretAccessKey']), len(c['sessionToken'])))

    hits = earthaccess.search_data(short_name='NISAR_L2_GCOV_PROVISIONAL_V1',
                                   granule_name=granule + '*')
    s3url = None
    for h in hits:
        for link in h.data_links(access='direct'):
            if link.startswith('s3://') and link.endswith('.h5'):
                s3url = link
    if s3url is None:
        print('no direct S3 link from CMR'); return 1
    bucket, key = s3url[5:].split('/', 1)
    print('object: s3://%s/%s' % (bucket, key))

    print('\n1. boto3 GetObject (1 byte) -- tests credentials and object permission')
    try:
        import boto3
        s3 = boto3.client('s3', region_name='us-west-2')
        r = s3.get_object(Bucket=bucket, Key=key, Range='bytes=0-0')
        print('   OK, read %d byte(s)' % len(r['Body'].read()))
        getOK = True
    except Exception as e:
        print('   FAILED: %s: %s' % (type(e).__name__, str(e)[:140]))
        getOK = False

    print('\n2. boto3 ListObjectsV2 -- GDAL may need this, ASF often does not grant it')
    try:
        boto3.client('s3', region_name='us-west-2').list_objects_v2(
            Bucket=bucket, Prefix=key, MaxKeys=1)
        print('   OK (ListBucket permitted)')
    except Exception as e:
        print('   DENIED: %s -- expected; GDAL must not rely on listing' % type(e).__name__)

    print('\n3. GDAL /vsis3, with settings that avoid listing')
    from osgeo import gdal
    gdal.UseExceptions()
    for name, opts in (
            ('as-is', {}),
            ('DISABLE_READDIR=YES', {'GDAL_DISABLE_READDIR_ON_OPEN': 'YES'}),
            ('+ virtual hosting', {'GDAL_DISABLE_READDIR_ON_OPEN': 'YES',
                                   'AWS_VIRTUAL_HOSTING': 'TRUE'}),
            ('+ regional endpoint', {'GDAL_DISABLE_READDIR_ON_OPEN': 'YES',
                                     'AWS_S3_ENDPOINT': 's3.us-west-2.amazonaws.com'}),
    ):
        for k, v in opts.items():
            gdal.SetConfigOption(k, v)
        try:
            st = gdal.VSIStatL('/vsis3/%s/%s' % (bucket, key))
            print('   %-22s stat OK, size %s' % (name, st.size if st else None))
        except Exception as e:
            print('   %-22s FAILED: %s' % (name, str(e)[:90]))
        for k in opts:
            gdal.SetConfigOption(k, None)
    if getOK:
        print('\nboto3 could read the object, so credentials are fine; the failures above are '
              'GDAL configuration.\nWhichever line says "stat OK" names the setting to use.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
