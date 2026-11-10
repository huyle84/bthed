# Stage 2a of BTHED. Requires the pq-crystals ML-DSA reference code:
#   git clone https://github.com/pq-crystals/dilithium.git ../../dilithium   (or reuse ../../stage1/dilithium)
REF ?= $(firstword $(wildcard ../../dilithium/ref ../../stage1/dilithium/ref))
CC ?= gcc
CFLAGS ?= -O3 -march=native -fomit-frame-pointer -Wall -Wextra -Wno-unused-parameter
SRC = $(REF)/sign.c $(REF)/packing.c $(REF)/polyvec.c $(REF)/poly.c $(REF)/ntt.c \
      $(REF)/reduce.c $(REF)/rounding.c $(REF)/fips202.c $(REF)/symmetric-shake.c \
      $(REF)/randombytes.c

all: bthed_fhe2 bthed_fhe3 bthed_fhe5 bthed_cmp

bthed_cmp: bthed_cmp.c
	$(CC) $(CFLAGS) -o $@ bthed_cmp.c -lm

bthed_fhe2 bthed_fhe3 bthed_fhe5: bthed_fhe.c $(SRC)
	$(CC) $(CFLAGS) -I$(REF) -DDILITHIUM_MODE=$(subst bthed_fhe,,$@) -o $@ bthed_fhe.c $(SRC) -lm

clean:
	rm -f bthed_fhe2 bthed_fhe3 bthed_fhe5 bthed_cmp
