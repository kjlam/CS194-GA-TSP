NVCC = nvcc


tsp.o: tsp.cu
	$(NVCC) -lcuda -lcudart ./tsp.cu -o tsp

tsp: tsp.o
	

