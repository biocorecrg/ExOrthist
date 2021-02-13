// base.tf

variable "profile" {
  type = string
  default = default
}

variable "credentials" {
  type = string
}

variable "region" {
  type = string
}

variable "amientrypoint" {
  type = string
}

variable "key_name" {
  type = string
}

variable "instance_entry" {
  type = string
}

provider "aws" {
  profile = var.profile
  shared_credentials_file = var.credentials
  region     = var.region
}

// You may define an entry point for convenience

resource "aws_instance" "entrypoint" {

  ami         = var.amientrypoint
  instance_type = var.instance_entry
  iam_instance_profile = "S3access"
  key_name = var.key_name
  security_groups = [ "allow_ssh" ]
  tags = {
	 name = "entrypoint"
  }

}

resource "aws_s3_bucket" "data-nf" {
  bucket = "data-nf"
  acl    = "private"
  force_destroy = true

  tags = {
    name = "s3-data-nf"
  }
}
