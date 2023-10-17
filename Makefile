PROGNAME:=advection_solver

DATADIR:=data/
PLOTSDIR:=plots/
VIDEODIR:=video/

DATA:=$(DATADIR)*.bin
PLOTS:=$(PLOTSDIR)*.png
VIDEO:=$(VIDEODIR)*.mp4

CC:=gcc
CFLAGS+=-fopenmp
LDLIBS+=-lm -pthread

SRC:=src/$(PROGNAME).c src/utils.c

.PHONY: clean purge setup video

advection_solver: $(SRC)
	$(CC) $^ $(CFLAGS) $(LDLIBS) -o $@

clean:
	-rm -f advection_solver

purge:
	-rm -f advection_solver $(DATA) $(PLOTS) $(VIDEO)

setup:
	-mkdir -p data plots video

video:
	./plot.sh > /dev/null
	ffmpeg -y -i $(PLOTSDIR)%5d.png -vf format=yuv420p $(VIDEODIR)animation.mp4 &> /dev/null