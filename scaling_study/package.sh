#!/bin/sh
tar cf macsio_scaling_study_final.tar package.sh make_jobs.py \
    results.dat results.html graph*.png \
    `find macsio_jobs -name macsio-timings.log` \
    `find macsio_jobs -name macsio-log.log`

