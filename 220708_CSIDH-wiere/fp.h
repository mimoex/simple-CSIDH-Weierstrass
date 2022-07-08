#pragma once
#include <iostream>
#include <gmpxx.h>
#include <random>
#include <chrono>

using namespace std;

const size_t N = 74;	//log(11^74)/log(2) à 255.9979

const int primes[N] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
						67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131,
						137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
						199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271,
						277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353,
						359, 367, 373, 587 };

struct Point {
	mpz_class x, y;
	bool inf = false;
};

//—LŒÀ‘Ì‚Ì‰ÁZ	c=a+b%p
void add_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//—LŒÀ‘Ì‚ÌŒ¸Z	c=a-b%p
void sub_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//—LŒÀ‘Ì‚ÌæZ	c=a*b%p
void mul_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//—LŒÀ‘Ì‚Ì‚×‚«æ	c=a^b%p
void pow_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//—LŒÀ‘Ì‚Ì‹tŒ³	c=a^-1
void inv_fp(const mpz_class& a, const mpz_class& p, mpz_class* c);
//—LŒÀ‘Ì‚ÌœZ	c=a/b%p
void div_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//pˆÈ‰º‚Ì—”‚ğ¶¬
mpz_class random_fp(const mpz_class& p);


Point copy_point(const Point& p);

//Weierstrass‹Èü‚Ì‰ÁZŒö®
Point ec_add(const Point& p, const Point& q, const mpz_class& a, const mpz_class& mod);

//‘È‰~‹Èüã‚Ì“_‚ğ¶¬
Point gen_point(const mpz_class& a, const mpz_class& b, const mpz_class& mod);
size_t check_point(const Point& p, const mpz_class a, const mpz_class b, const mpz_class& mod);
Point gen_point_sqrt(const mpz_class& a, const mpz_class& b, const mpz_class& mod);
