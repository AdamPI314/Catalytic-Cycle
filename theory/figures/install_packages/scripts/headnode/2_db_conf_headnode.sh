#!/bin/sh

stty sane
echo /tmp/db_conf.log.$$ > /tmp/db_conf.log.$$ 2>&1

sudo update-rc.d mysql enable >> /tmp/db_conf.log.$$ 2>&1
sudo service mysql start >> /tmp/db_conf.log.$$ 2>&1

# do the following manually, basic database setup, create a dababase and a table, grant privileges
# notice the second line, the privileges will granted, if the user doesn't exist, it will be
# created first
sudo mysql -u root <<EOF
create database if not exists slurm_acct_db;
grant all privileges on slurm_acct_db.* to 'slurm'@'localhost';
set password for 'slurm'@'localhost' = password('slurmdbpass');
grant usage on *.* to 'slurm'@'localhost';
grant all privileges on slurm_acct_db.* to 'slurm'@'localhost';
flush privileges;
exit
EOF
