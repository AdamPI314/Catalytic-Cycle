#!/bin/sh

stty sane
echo /tmp/start_compute.log.$$ > /tmp/start_compute.log.$$ 2>&1

sudo systemctl enable munge >> /tmp/start_compute.log.$$ 2>&1
sudo systemctl start munge >> /tmp/start_compute.log.$$ 2>&1

sudo systemctl enable slurmd >> /tmp/start_compute.log.$$ 2>&1
sudo systemctl start slurmd >> /tmp/start_compute.log.$$ 2>&1
