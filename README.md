# data-portal-summary-stats

The service generates per-project summary figures. The service updates daily by running 
a Docker container on AWS Fargate. Figures are persisted on AWS S3. The general
architecture is as follows:

![](./illustrations/spec_v4.png)

## Set-up
The code is written for Python version 3.6. Clone the repository to your local system and 
navigate into the `data-portal-summary-stats` directory. Create a virtual environment
and run `pip install -r requirements.txt` to install dependencies. 

#### Install AWS
The service runs on AWS. Install the [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/install-linux.html)
 and have your credentials and configuration files in `~/.aws`. Authentication and access is role-based. 

#### Installation and configuration of Terraform
This app uses the infrastructure management software [Terraform](https://learn.hashicorp.com/terraform/getting-started/install.html).
 Download the binary executable (version 0.12 
or higher) into a directory that is part of the system's Linux `$PATH` so you can execute it 
anywhere. Run `terraform init` to create the Terraform backend. In the directory `infra` 
create a file `variables.tf` which contains all required variables. This table contains
a list of required variable names and their recommended values (all values are of type _String_).
 

| Variable | Description | Value |
| --- | --- | --- |
| `dpss_task_memory` | RAM allocated to the container instance [GB] | "16384" |
| `dpss_task_cpu` | The number of CPU units (1024 is equivalent to 1 core)| "2048" |



### Build the Docker image
To build the Docker image runs `data-portal-summary-stats` create a directory that you
name for instance `build`, then copy only the relevant files from the project root into `build`,
such that the directory tree looks like so:
```bash
.
├── data_portal_summary_stats.py
├── Dockerfile
├── requirements.txt
└── src
    ├── matrix_summary_stats.py
    ├── runner.py
    ├── settings.py
    └── utils.py

1 directory, 7 files
```


Docker commands need to be run in sudo mode for security reasons. You can avoid this by adding
 the user to the docker group, see [here](https://linoxide.com/linux-how-to/use-docker-without-sudo-ubuntu/))
 Next navigate to the `build` directory, prepared as explained in the previous paragraph, and
 execute:
```bash
docker build --tag=data-portal-summary-stats:<tag #> .
```
where `tag #` denotes a dot version of the build (e.g., `0.1`).
#### Running the `data-portal-summary-stats` Docker container

This describes how to run the images as a container on your local system. The code in 
`data_portal_summary_stats` needs to access AWS resources, so you need to pass your AWS
 credentials to container run. Suppose your credentials are in the usual subdirectory `~/.aws` 
 in your home directory, and they contain a profile named `my-profile`. In that case mount that 
 subdirectory as a volume inside the container, and set the environment variable with the default 
 profile like so:
```bash
docker run -v /home/michael/.aws:/root/.aws -e AWS_DEFAULT_PROFILE=my-profile \
       data-portal-summary-stats:<tag #> --environ dev --source fresh --blacklist true
```
The three arguments, `--environ`, `--source`, and `--blacklist` and their values follow the name 
of the image and its tag. 


### Deploy the service
To deploy with Terraform navigate to `infra` from the project root, then run 
```bash
export AWS_PROFILE=your-profile
```
followed by `terraform plan`, followed by `terraform apply` to deploy the infrastructure. This sets
up a `Scheduled Task`

### Handling of large matrix files
Depending on some parameters the matrix files of some projects are too large to process for even the largest available ECS 
container instance (e.g., we found that a matrix file of > 2 GB cannot be processed). We offer
 two solutions to this problem. 
 
The **first** and preferred solution aims to
 process all matrices by choosing a combination of RAM and a setting a suitable 
 value for the input
 argument `--min_gene_count`. The following figure shows results of some tests we ran using
 a matrix of project _Census of Immune Cells_ on September 4, 2019. At this point this is
 the largest matrix files of the HCA. The blue area denotes processing failures, and the 
 orange area denotes processing successes:
  
  ![tests](./illustrations/large-file-experiment_figure.png)
  
Based on these results we chose a combination of 16 GB RAM and setting `min_gene_count` to 
a value of `1200` as the optimal combination. 
 
The **second** solution is to simply black-list the matrix files that are too large to
 process. _blacklisting_
 a matrix file excludes it from the processing and is a solution of last resort. To use 
 the blacklist feature, set the command input argument `--blacklist` to `true` and create a 
 file named `blacklist`. To specify a matrix file from processing add its project ID to 
 that file, following the style of one project ID per line, without any delimiters. Next
 upload the file `blacklist` to the root of the S3 bucket of the corresponding
 deployment environment. One way to upload the file is (assuming `blacklist` is in the 
 directory :
 ```bash
 ls blacklist   // should return `blacklist`
 export AWS_DEFAULT_PROFILE=my_profile
 aws s3 cp blacklist s3://my_bucket/ 
```
where _my_bucket_ is the bucket name of the deployment environment.  