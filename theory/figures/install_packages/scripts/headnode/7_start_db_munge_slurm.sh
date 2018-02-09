#!/bin/sh

stty sane
echo /tmp/start_head.log.$$ > /tmp/start_head.log.$$ 2>&1

sudo update-rc.d mysql enable >> /tmp/start_head.log.$$ 2>&1
sudo service mysql start >> /tmp/start_head.log.$$ 2>&1

sudo systemctl enable munge >> /tmp/start_head.log.$$ 2>&1
sudo systemctl start munge >> /tmp/start_head.log.$$ 2>&1

sudo systemctl daemon-reload >> /tmp/start_head.log.$$ 2>&1 
sudo systemctl enable slurmdbd >> /tmp/start_head.log.$$ 2>&1
sudo systemctl start slurmdbd >> /tmp/start_head.log.$$ 2>&1
sudo systemctl enable slurmctld >> /tmp/start_head.log.$$ 2>&1
sudo systemctl start slurmctld >> /tmp/start_head.log.$$ 2>&1
