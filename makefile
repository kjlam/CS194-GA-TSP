NVCC = nvcc

tsp: tsp.o
	$(NVCC) -lcuda -lcudart ./tsp.cu -o tsp 


tsp.o: tsp.cu
	$(NVCC) -lcuda -lcudart ./tsp.cu -o tsp
