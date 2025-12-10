nohup taskset -c 0-3 bash rfd_binder_1.sh &
nohup taskset -c 4-7 bash rfd_binder_2.sh &
nohup taskset -c 8-11 bash rfd_binder_3.sh &
nohup taskset -c 12-15 bash rfd_binder_4.sh &
