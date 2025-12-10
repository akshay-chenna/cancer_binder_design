<< 'END'
for i in {1..9}
do
	bash rfd_binder_partial_${i}.sh &
	bash rfd_binder_partial_$((i+8)).sh &
	wait
done 
END

for i in {6..8}
do
	bash rfd_binder_partial_${i}.sh &
done 

for i in {14..16}
do
	bash rfd_binder_partial_${i}.sh &
done

for i in {17..18}
do
	bash rfd_binder_partial_${i}.sh &
done

wait
