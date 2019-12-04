class Config:
    def __init__(self, stage):
        assert stage in ('dev', 'integration', 'staging', 'prod')
        self.stage = stage

    @property
    def stage_str(self):
        return '' if self.stage == 'prod' else f'{self.stage}.'

    @property
    def azul_project_endpoint(self):
        return f'https://service.{self.stage_str}explore.data.humancellatlas.org/repository/projects/'

    @property
    def hca_matrix_service_endpoint(self):
        return f'https://matrix.{self.stage_str}data.humancellatlas.org/v1/'

    @property
    def s3_bucket_name(self):
        return f'{self.stage_str}project-assets.data.humancellatlas.org'

    @property
    def s3_canned_matrix_prefix(self):
        return 'project-assets/project-matrices/'

    @property
    def s3_figures_prefix(self):
        return 'project-assets/project-stats/'
