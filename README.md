# data-portal-summary-stats

The application generates the following per-project summary figures:

1. List of genes with highest cell count
2. List of top 20 genes with the highest expression variability
 
The application updates daily by running a Docker container on AWS Fargate. Figures are persisted 
on AWS S3. The general architecture is as follows:

![](./illustrations/spec_v4.png)

## Set-up
_data-portal-summary-stats_ is written in Python version 3.6. Clone the repository to your 
local system and navigate into the `data-portal-summary-stats` directory. If you wish to run code
from source create a virtual environment, name it `.venv`, and run 
`pip install -r requirements.txt` to install dependencies. 
But the intention is to run the code in a Docker container. So 
[Docker](https://www.docker.com) needs to be installed on your system. The application 
uses AWS. Install AWS CLI and have the credentials configured for AWS. 

#### Installation and configuration of Terraform
This app uses the infrastructure management software 
[Terraform](https://learn.hashicorp.com/terraform/getting-started/install.html), version 0.12 
or higher. In the directory `infra` run `terraform init` to create the Terraform backend. 
Next source the environment file for the specific deployment environment you would like to set 
up (e.g., _dev_). From the project root run 
(requires the utility [`jq`](https://stedolan.github.io/jq/)):

```bash
source environment.dev
```

## Running the application

Once all required software is installed, the basic steps to start the application are
1. Build the Docker image
2. Push the image to the AWS registry
3. Deploy the application using the Terraform

In the following we give detailed instructions for each individual step.

### 1. Build the Docker image
To build the Docker image `data-portal-summary-stats` execute

```bash
docker build --tag=data-portal-summary-stats:$TF_VAR_image_tag .
```
from the project root.

#### Input arguments to the container

Once the image is built on your local system run 
`docker run data-portal-summary-stats:$I$TF_VAR_image_tag -h` to 
return the help message with a brief description of the command line arguments. For more 
information on the last argument, `--min_gene_count`, 
[go the Matrix Service Swagger UI](https://matrix.staging.data.humancellatlas.org/)
and exercise the endpoint `/v1/filters/{filter_name}` with `genes_detected` as `filter_name`.   

#### Running the `data-portal-summary-stats` Docker container locally

To run the image as a container on your local system pass your AWS
 credentials to it. Suppose your credentials are in the usual subdirectory `~/.aws` 
 in your home directory, and they contain a profile named `my-profile`. In that case mount that 
 subdirectory as a volume inside the container, and set the environment variable with the default 
 profile like so:
 
```bash
docker run -v $HOME/.aws:/root/.aws -e AWS_DEFAULT_PROFILE=$AWS_DEFAULT_PROFILE \
       data-portal-summary-stats:$TF_VAR_image_tag \ 
       --environ $TF_VAR_dev --source $TF_VAR_source --min_gene_count $TF_VAR_min_gene_count
```

### 2. Push the image to AWS ECR
First [create a repository](https://console.aws.amazon.com/ecr/repositories) for 
your images using the AWS console. Next push the image you just created to 
[ECR](https://aws.amazon.com/ecr/), AWS's image registry. Be sure to give that repository the 
same name as the image. Authentication is needed to push any image from localhost to that registry. 
Open a terminal and run

```bash
aws ecr get-login --region $TF_VAR_aws_region --no-include-email
```
This prints the Docker log-in command, something like 
`docker login -u AWS -p password https://$TF_VAR_acc_number.dkr.ecr.$TF_VAR_aws_region.amazonaws.com`, where 
`password` has several hundred characters. Copy this command and run it in the terminal. Then 
tag the image using the same tag in the repository by running

```bash
docker tag data-portal-summary-stats:$TF_VAR_image_tag $TF_VAR_acc_number.dkr.ecr.$TF_VAR_aws_region.amazonaws.com/data-portal-summary-stats:$TF_VAR_image_tag
```
and push it to the repository to the created namespace:

```bash
docker push $TF_VAR_acc_number.dkr.ecr.$TF_VAR_aws_region.amazonaws.com/data-portal-summary-stats:$TF_VAR_image_tag
```

### 3. Deploy the application
We use AWS's ECS container orchestration service with the Fargate launchtype to run the 
container as a _task_ and Terraform to deploy and manage the application. Next, navigate to 
`infra` and run `terraform plan` if you want to inspect what resources will be created. Otherwise, 
deploy the infrastructure and start the application by running `terraform apply`. Follow the 
instructions on the screen. It should end on "Apply complete! Resources [...]". If you see that 
message the deployment process is complete. Once deployed, the application runs on a schedule to 
create a set of summary statistics figures every 24 hours and can be monitored on AWS CloudWatch.

## Handling of large matrix files
Depending on some parameters the matrix files of some projects are too large to process for even 
the largest available ECS 
container instance (e.g., we found that a matrix file of > 2 GB size cannot be processed). We offer
 two solutions to this problem. 
 
The **first** and preferred solution aims to
 process all matrices by choosing a combination of sufficient RAM and setting 
  input argument `--min_gene_count` to a suitable value. The large that value, the smaller the
  matrix file size. The following figure shows results of some tests we ran using
 a matrix of project _Census of Immune Cells_ from September 4, 2019, which was the project with 
 the largest matrix file in the HCA on that day. The orange area denotes successful processing 
 and the blue area denotes failed processing:
  
  ![tests](./illustrations/large-file-experiment_figure.png)
  
Based on these results we chose a combination of 16 GB RAM and setting `min_gene_count` to 
a value of `1200` as the optimal combination. That creates matrix files of a size such that
all can be processed.
 
The **second** solution is to simply black-list the matrix files that are too large to
 process. _blacklisting_
 a matrix file excludes it from the processing and is a solution of last resort. To use 
 the blacklist feature, set the command input argument `--blacklist` and create a 
 file named `blacklist`. To specify a matrix file from processing add its project ID to 
 that file, following the style of one project ID per line, without any delimiters. Next
 upload the file `blacklist` to the root of the `project-assets` S3 bucket of the corresponding
 deployment environment.
 
 ## Tests
 Run `python -m unittest` from the project root.
 