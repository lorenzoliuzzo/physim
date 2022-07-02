compile all: harmonic_oscill.x forced_oscill.x damped_oscill.x forced_damped_oscill.x

esegui all: compile all
	./harmonic_oscill.x
	./forced_oscill.x
	./damped_oscill.x
	./forced_damped_oscill

harmonic_oscill.x: harmonic_oscill.cpp
	g++ harmonic_oscill.cpp -o harmonic_oscill.x


forced_oscill.x: forced_oscill.cpp
	g++ forced_oscill.cpp -o forced_oscill.x


damped_oscill.x: damped_oscill.cpp
	g++ damped_oscill.cpp -o damped_oscill.x


forced_damped_oscill.x: forced_damped_oscill.cpp
	g++ forced_damped_oscill.cpp -o forced_damped_oscill.x

clean:
	rm *.x