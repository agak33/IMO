
compile:
	g++ main.cpp -o out

run:
	mkdir -p output plots
	./out
	cd ..; python visualisation.py -i local_search/output -f greedy_tsp/input -o local_search/plots

viz:
	cd ..; python visualisation.py -i local_search/output -f greedy_tsp/input -o local_search/plots

clean:
	rm ./output/*
	rm ./plots/*