rm results

./bin/roadhog --alg ch --problem ./dimacs/inputs/USA-road-d/USA-road-d.NY.p2p --input ./dimacs/USA-road-d.NY.gr.ch ./dimacs/USA-road-d.NY.co ./dimacs/USA-road-d.NY.gr.ooc | tee -a results

./bin/roadhog --alg ch-astar --problem ./dimacs/inputs/USA-road-d/USA-road-d.NY.p2p --input ./dimacs/USA-road-d.NY.gr.ch ./dimacs/USA-road-d.NY.co ./dimacs/USA-road-d.NY.gr.ooc | tee -a results

./bin/roadhog --alg fch-bbaf --problem ./dimacs/inputs/USA-road-d/USA-road-d.NY.p2p --input ./dimacs/USA-road-d.NY.gr.ch ./dimacs/USA-road-d.NY.co ./dimacs/USA-road-d.NY.gr.ooc ./dimacs/USA-road-d.NY.gr.ch.fch-bbaf.arclabel ./dimacs/USA-road-d.NY.gr.metis.part.128 | tee -a results

./bin/roadhog --alg fch-bb --problem ./dimacs/inputs/USA-road-d/USA-road-d.NY.p2p --input ./dimacs/USA-road-d.NY.gr.ch ./dimacs/USA-road-d.NY.co ./dimacs/USA-road-d.NY.gr.ooc ./dimacs/USA-road-d.NY.gr.ch.fch-bb.arclabel ./dimacs/USA-road-d.NY.gr.metis.part.128 | tee -a results

./bin/roadhog --alg fch-af --problem ./dimacs/inputs/USA-road-d/USA-road-d.NY.p2p --input ./dimacs/USA-road-d.NY.gr.ch ./dimacs/USA-road-d.NY.co ./dimacs/USA-road-d.NY.gr.ooc ./dimacs/USA-road-d.NY.gr.ch.fch-af.arclabel ./dimacs/USA-road-d.NY.gr.metis.part.128 | tee -a results

./bin/roadhog --alg fchx --problem ./dimacs/inputs/USA-road-d/USA-road-d.NY.p2p --input ./dimacs/USA-road-d.NY.gr.ch ./dimacs/USA-road-d.NY.co ./dimacs/USA-road-d.NY.gr.ooc | tee -a results

./bin/roadhog --alg astar --problem ./dimacs/inputs/USA-road-d/USA-road-d.NY.p2p --input ./dimacs/USA-road-d.NY.gr ./dimacs/USA-road-d.NY.co ./dimacs/USA-road-d.NY.gr.ooc | tee -a results

./bin/roadhog --alg bi-astar --problem ./dimacs/inputs/USA-road-d/USA-road-d.NY.p2p --input ./dimacs/USA-road-d.NY.gr ./dimacs/USA-road-d.NY.co ./dimacs/USA-road-d.NY.gr.ooc | tee -a results

./bin/roadhog --alg dijkstra --problem ./dimacs/inputs/USA-road-d/USA-road-d.NY.p2p --input ./dimacs/USA-road-d.NY.gr ./dimacs/USA-road-d.NY.co ./dimacs/USA-road-d.NY.gr.ooc | tee -a results

./bin/roadhog --alg bi-dijkstra --problem ./dimacs/inputs/USA-road-d/USA-road-d.NY.p2p --input ./dimacs/USA-road-d.NY.gr ./dimacs/USA-road-d.NY.co ./dimacs/USA-road-d.NY.gr.ooc | tee -a results


