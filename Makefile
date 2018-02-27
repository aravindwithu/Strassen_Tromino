run:
	gcc strassen.c -lm -o strassen.out
	gcc tromino.c -lm -o tromino.out

clean:
	rm strassen.out
	rm tromino.out
