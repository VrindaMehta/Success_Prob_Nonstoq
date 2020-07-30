CC = icc
#CFLAGS= -pedantic -Wall -Weffc++ -Wextra -I/usr/include/openblas -lopenblas

succprob: Essential.cpp Read.cpp Parameters_Nonstoq.cpp Product_Nonstoq.cpp Lanczos_succ.cpp QA.cpp
	$(CC) -qopenmp -O3 Essential.cpp Read.cpp Parameters_Nonstoq.cpp Lanczos_succ.cpp QA.cpp Product_Nonstoq.cpp -o lanczos_prob -mkl

instover: Essential.cpp Read.cpp Parameters_Nonstoq.cpp Product_Nonstoq.cpp StateOverlap.cpp QA_StateOverlap.cpp 
	$(CC) -O3 Essential.cpp Read.cpp Parameters_Nonstoq.cpp Product_Nonstoq.cpp StateOverlap.cpp QA_StateOverlap.cpp -o lanczos_instoverlap -mkl
clean: 
	rm -f *.o lanczos_prob lanczos_instoverlap 
