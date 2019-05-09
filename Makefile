object_MC = main.o tools.o initialize.o

all: mc

mc: $(object_MC)
	gcc -o mc $(object_MC) -lm -g

main.o: main.c tools.h
	gcc -c main.c -g

tools.o: tools.c tools.h
	gcc -c tools.c -g

initialize.o: initialize.c initialize.h
	gcc -c initialize.c -g	

clean:
	rm *.o	
