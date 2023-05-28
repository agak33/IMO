
compile:
	g++ main.cpp -o out

run:
	mkdir -p output plots
	./out
	cd ..; python visualisation.py -i genetic/output -f genetic/input -o genetic/plots

viz:
	cd ..; python visualisation.py -i genetic/output -f genetic/input -o genetic/plots

clean:
	rm ./output/*
	rm ./plots/*