all:
	g++ -o generate_field circu_embed.cpp -lfftw3 -lm covariance_fctn.cpp
