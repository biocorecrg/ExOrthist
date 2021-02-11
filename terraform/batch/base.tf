// Credentials to use. Here in Frankfurt region
// Different ways to store credentials
// provider "aws" {
//  access_key = "yourkeyhere"
//  secret_key = "yoursecrethere"
//  region     = "eu-central-1"
// }

provider "aws" {
  profile = "default"
  shared_credentials_file = "/home/myusername/.aws/credentials"
  region     = "eu-central-1"
}

// You may define an entry point for convenience

resource "aws_instance" "entrypoint" {

  ami         = "ami-06a3c664bd6bb3fce"
  instance_type = "t2.micro"
  iam_instance_profile = "S3access"
  key_name = "key-nf"
  security_groups = [ "allow_ssh" ]
  tags = {
	name = "entrypoint"
  }

}

resource "aws_s3_bucket" "exorthist-nf" {
  bucket = "exorthist-nf"
  acl    = "private"
  force_destroy = true

  tags = {
    name = "S3 for ExOrthist"
  }
}
