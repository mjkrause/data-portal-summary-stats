# data-portal-summary-stats

The application generates the following per-project summary figures:

1. Highest expressing genes
2. Violin plots of cells, all genes, and percent of mitochondrial genes
3. Number of genes over number of counts.
4. Percent mitochondrial genes over number of counts.
5. Visualize highly-variable genes
6. Principal components, PC2 against PC1
7. tSNE, Umap 2 against Umap1, of Louvain and CST3.
8. Ranks genes groups.

The application updates daily by running a Docker container on AWS Fargate. Figures are persisted 
on AWS S3. Figures will be assigned keys based on the matrix project ID and the 
library construction approach. Matrices containing cells processed using both 10X sequencing
and Smart-seq2 will be split into separate entities before analysis.

An example of a figure key is 

`project-assets/project-stats/abe1a013-af7a-45ed-8c26-f3793c24a1f4.homo_sapiens/10X/violin.png`

For a canned matrix called `abe1a013-af7a-45ed-8c26-f3793c24a1f4.homo_sapiens.mtx.zip`

The general architecture is as follows:

![](./illustrations/spec_v4.png)

## Set-up
_data-portal-summary-stats_ is written in Python version 3.6. Clone the repository to your 
local system and navigate into the `data-portal-summary-stats` directory. 
If using PyCharm, mark the `src` and `test` directories as source roots.
If you wish to run code
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
source environment.sh
```

## Running the application

Once all required software is installed, the basic steps to start the application are
1. Build the Docker image
2. Push the image to the AWS registry
3. Deploy the application using the Terraform

In the following we give detailed instructions for each individual step.

#### Input arguments to the container

The project includes two files responsible for setting environment variables 
used to configure the application. `environment.env` defines the variables 
needed at runtime in a format suitable for the EnvFile Pycharm extension.
`environment.sh` exports the variables defined in `environment.env` as well as 
additional variables needed by docker and terraform during building and 
deployment.

`environemnt.sh` should be `source`-ed before building the container.
The variables are documented in `environemnt.env`. those most important for configuration are:

 - `DPSS_MTX_SOURCE_STAGE`: Which deployment to acquire resources from, such as 
 blacklists and matrices. Should be set to one of `dev`, `integration`, `staging`, or `prod`.

 - `DPSS_MTX_TARGET_STAGE`: Which deployment bucket to upload generated figures to.
 Should be one of `dev`, `integration`, `staging`, `prod`, or `ux-dev`.

 - `DPSS_MATRIX_SOURCE`: Whether to generate stats for "fresh" (retrieved from 
 HCA matrix service) or "canned" (stored in S3) matrices. Should be either `fresh` or `canned`.

 - `DPSS_BLACKLIST`: Whether or not to use the blacklist file to filter matrices.
 Set to `1` to use blacklist; another other value or left unset to ignore it.

### 1. Build the Docker image

To build the Docker image `data-portal-summary-stats` execute

```bash
docker build --tag=data-portal-summary-stats:$DPSS_IMAGE_TAG .
```
from the project root.


#### Running the `data-portal-summary-stats` Docker container locally

To run the image as a container on your local system pass your AWS
 credentials to it. Suppose your credentials are in the usual subdirectory `~/.aws` 
 in your home directory, and they contain a profile named `my-profile`. In that case mount that 
 subdirectory as a volume inside the container, and set the environment variable with the default 
 profile like so:
 
```bash
docker run -v $HOME/.aws:/root/.aws --env-file environment.env -e AWS_DEFAULT_PROFILE=$AWS_DEFAULT_PROFILE data-portal-summary-stats:$TF_VAR_image_tag
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
 
The current implementation attempts to process all matrices by choosing a 
combination of sufficient RAM and dropping rows from the matrices based on the
number of genes detected. However, this is not a viable long-term solution
as the gene threshold parameter affects the analysis.

The following figure shows results of some tests we ran using a matrix of project 
_Census of Immune Cells_ from September 4, 2019, which was the project with 
 the largest matrix file in the HCA on that day. The orange area denotes 
 successful processing and the blue area denotes failed processing:
  
  ![tests](./illustrations/large-file-experiment_figure.png)
  
Based on these results we chose a combination of 16 GB RAM and setting `min_gene_count` to 
a value of `1200` as the optimal combination. That creates matrix files of a size such that
all can be processed.

1200 is currently still the threshold used for all cells and matrices regardless 
of library construction approach(es). However, we have not fully investigated
whether this value is truly appropriate for the analysis. As such, other measures
to address large matrices may be needed in the future as these thresholds are tweaked.  
 
The second solution is to simply blacklist the matrix files that are too large to
 process. _blacklisting_
 a matrix file excludes it from the processing and is a solution of last resort. To use 
 the blacklist feature, set the environment variable `blacklist=1` and create a 
 file named `blacklist`. To specify a matrix file from processing add its project ID to 
 that file, following the style of one project ID per line, without any delimiters. Next
 upload the file `blacklist` to the root of the `project-assets` S3 bucket of the corresponding
 deployment environment.
 
 ## Tests
 Source `environment.sh` and run `python -m unittest` from the project root.
 If using PyCharm, consider installing the EnvFile extension and using it with `environemnt.env`.
 