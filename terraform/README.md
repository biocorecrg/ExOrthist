Example [Terraform](https://www.terraform.io/) recipe for running an AWS Batch with CPU node instances

* Entrypoint AMI can be created following instructions here: https://www.nextflow.io/docs/latest/awscloud.html#custom-ami

Terraform variables values can be stored as exported BASH variables, such as:

```
export TF_VAR_amientrypoint=ami-xxxxxx
export TF_VAR_amibatch=ami-xxxxx
export TF_VAR_key_name=mykeyname
export TF_VAR_instance_entry=t2.micro
export TF_VAR_region=eu-central-1
export TF_VAR_profile=default
export TF_VAR_credentials=/home/myuser/.aws/credentials

```

or otherwise you will be prompted when running Terraform.

A good practice is keeping them in a file (e.g. <code>myvars.env</code>) and declare them by ```source myvars.env``` before running the commands below.

## Instructions

Go to the directory with the customized Terraform scripts and check they are valid:

```terraform validate```

You may need to ```terraform init``` first.

To review which infrastructure will be deployed:

```terraform plan```

To actually deploy it:

```terraform apply```

Then you need to connect to entrypoint to execute Nextflow. A convenient example SSM command is also provided under ```run.sh``` if that method is preferred.

Once work is done:

```terraform destroy```

* *Be careful, configuration is set up to destroy S3 Bucket by default and you might want to change that.*

### Uploading content to S3 Backet

It is very straightforward to upload contents to your S3 bucket.

An example command:

    aws s3 cp --recursive . s3://data-nf/exorthist/input/benchmark

You don't need to create target directories in advance.

## Additional links
* [Storing sensitive data with Terraform](https://blog.gruntwork.io/a-comprehensive-guide-to-managing-secrets-in-your-terraform-code-1d586955ace1)
* [AWS CLI commands for S3 buckets](https://docs.aws.amazon.com/cli/latest/userguide/cli-services-s3-commands.html)
