
compile:
	g++ main.cpp -o out

run:
	mkdir -p output plots
	./out
	cd ..; python visualisation.py -i local_extension/output -f local_extension/input -o local_extension/plots

viz:
	cd ..; python visualisation.py -i local_extension/output -f local_extension/input -o local_extension/plots

clean:
	rm ./output/*
	rm ./plots/*