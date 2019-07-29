import boto3
import datetime
import logging
from chalice import Chalice, Rate
from chalicelib.settings import s3_bucket
from chalicelib.matrix_summary_stats import MatrixSummaryStats

app = Chalice(app_name='portal-summary-stats-lambda')


# _MSS = None  # matrix summary stats instance global
# _MSS_CLIENT = None
#
#
# def mss_client():
#     global _MSS_CLIENT
#     if _MSS_CLIENT is None:
#         _MSS_CLIENT = boto3.client('s3')
#     return _MSS_CLIENT
#
#
# def get_mss():
#     global _MSS
#     if _MSS is None:
#         _MSS = MatrixSummaryStats(s3_bucket['bucket_name'], s3_bucket['key'], mss_client)
#     return _MSS


@app.route('/')
def index():
    return {'sanity': 'OK'}


@app.schedule(Rate(1, unit=Rate.MINUTES))
def get_data_and_upload_images(event):

    #app.log.info('Downloading matrix, ETF, uploading images to S3...')
    # get_mss().get_expression_matrix()
    # get_mss().unzip_files()
    # get_mss().create_images()
    # get_mss().upload_figs_to_s3()

    bucket_name = s3_bucket['bucket_name']
    key = s3_bucket['key']

    t = str(datetime.datetime.now())
    #app.log.info(f'Uploading file at {t}')
    s3 = boto3.resource('s3')
    s3.Object(bucket_name, key + 'newfile.txt').put(Body=t)
