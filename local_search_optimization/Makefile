
compile:
	g++ main.cpp -o out

run:
	mkdir -p output plots
	./out
	cd ..; python visualisation.py -i local_search_optimization/output -f local_search_optimization/input -o local_search_optimization/plots

viz:
	cd ..; python visualisation.py -i local_search_optimization/output -f local_search_optimization/input -o local_search_optimization/plots

clean:
	rm ./output/*
	rm ./plots/*