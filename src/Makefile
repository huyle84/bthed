# Stage-1 BTHED simulator. Requires the pq-crystals ML-DSA reference code:
#   git clone https://github.com/pq-crystals/dilithium.git ../dilithium
REF ?= ../dilithium/ref
CC ?= gcc
CFLAGS ?= -O3 -march=native -fomit-frame-pointer -Wall -Wextra -Wno-unused-parameter
SRC = $(REF)/sign.c $(REF)/packing.c $(REF)/polyvec.c $(REF)/poly.c $(REF)/ntt.c \
      $(REF)/reduce.c $(REF)/rounding.c $(REF)/fips202.c $(REF)/symmetric-shake.c \
      $(REF)/randombytes.c

all: bthed_sim2 bthed_sim3 bthed_sim5

bthed_sim2 bthed_sim3 bthed_sim5: bthed_sim.c $(SRC)
	$(CC) $(CFLAGS) -I$(REF) -DDILITHIUM_MODE=$(subst bthed_sim,,$@) -o $@ bthed_sim.c $(SRC)

clean:
	rm -f bthed_sim2 bthed_sim3 bthed_sim5
