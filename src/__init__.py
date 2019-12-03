class Config:
    def __init__(self, environ):
        assert environ in ('dev', 'integration', 'staging', 'prod')
        self.environ = environ

    @property
    def environ_str(self):
        return '' if self.environ == 'prod' else f'{self.environ}.'

    @property
    def azul_project_endpoint(self):
        return f'https://service.{self.environ_str}explore.data.humancellatlas.org/repository/projects/'

    @property
    def hca_matrix_service_endpoint(self):
        return f'https://matrix.{self.environ_str}data.humancellatlas.org/v1/'

    @property
    def s3_bucket_name(self):
        return f'{self.environ_str}project-assets.data.humancellatlas.org'

    @property
    def s3_canned_matrix_prefix(self):
        return 'project-assets/project-matrices/'

    @property
    def s3_figures_prefix(self):
        return 'project-assets/project-stats/'
