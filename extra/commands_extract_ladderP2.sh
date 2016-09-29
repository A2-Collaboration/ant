#!/bin/bash

mysql -B -h slowcontrol.online.a2.kph -u report -p'$report' -D archive_backup\
    -e 'select unix_timestamp(smpl_time) ,  avg(float_val) FROM Run_2014_07_EPT_Prod WHERE channel_id=1870 GROUP BY Date(smpl_time),hour(smpl_time),minute(smpl_time);' \
    > data_jul

mysql -B -h slowcontrol.online.a2.kph -u report -p'$report' -D archive_backup\
    -e 'select unix_timestamp(smpl_time) ,  avg(float_val) FROM Run_2014_10_EPT_Prod WHERE channel_id=1870 GROUP BY Date(smpl_time),hour(smpl_time),minute(smpl_time);' \
    > data_okt

mysql -B -h slowcontrol.online.a2.kph -u report -p'$report' -D archive_backup\
    -e 'select unix_timestamp(smpl_time) ,  avg(float_val) FROM Run_2014_12_EPT_Prod WHERE channel_id=1870 GROUP BY Date(smpl_time),hour(smpl_time),minute(smpl_time);' \
    > data_dec
