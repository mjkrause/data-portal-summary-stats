endpoints = {
    'dev': {
        'hca_matrix_service_url': 'https://matrix.dev.data.humancellatlas.org/v1/',
        'azul_projects_url': 'https://service.dev.explore.data.humancellatlas.org/repository/projects/'
    },
    'integration': {
        'hca_matrix_service_url': 'https://matrix.integration.data.humancellatlas.org/v1/',
        'azul_projects_url': 'https://service.integration.explore.data.humancellatlas.org/repository/projects/'
    },
    'staging': {
        'hca_matrix_service_url': 'https://matrix.staging.data.humancellatlas.org/v1/',
        'azul_projects_url': 'https://service.staging.explore.data.humancellatlas.org/repository/projects/'
    },
    'prod': {
        'hca_matrix_service_url': 'https://matrix.data.humancellatlas.org/v1/',
        'azul_projects_url': 'https://service.explore.data.humancellatlas.org/repository/projects/'
    }
}

s3_bucket_info = {
    'dev': {
        'bucket_name': 'dev.project-assets.data.humancellatlas.org',
        'key': 'project-assets/project-stats/'
    },
    'integration': {
        'bucket_name': 'integration.project-assets.data.humancellatlas.org',
        'key': 'project-assets/project-stats/'
    },
    'staging': {
        'bucket_name': 'staging.project-assets.data.humancellatlas.org',
        'key': 'project-assets/project-stats/'
    },
    'prod': {
        'bucket_name': 'project-assets.data.humancellatlas.org',
        'key': 'project-assets/project-stats/'
    }

}

