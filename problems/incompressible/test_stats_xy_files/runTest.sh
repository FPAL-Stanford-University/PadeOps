#!/bin/bash

setup_padops
export CODEDIR=/home1/06632/ryanhass/codes/PadeOps/build/problems/incompressible

ibrun $CODEDIR/test_stats_xy input_interscale.dat blah | tee test_interscale.log
ibrun $CODEDIR/test_stats_xy input.dat | tee test64.log

sed -i "4,6s/64/128/g" input.dat
ibrun $CODEDIR/test_stats_xy input.dat | tee test128.log

sed -i "4,6s/128/256/g" input.dat
ibrun $CODEDIR/test_stats_xy input.dat | tee test256.log

sed -i "4,6s/256/64/g" input.dat

module load python
python gather_and_plot_error.py | tee process_error.log
