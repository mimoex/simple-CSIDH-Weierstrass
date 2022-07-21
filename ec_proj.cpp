#include "fp.h"

//https://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#doubling-dbl-2007-bl
Point_proj ec_double_proj(const Point_proj& p, const mpz_class& a, const mpz_class& mod)
{
	mpz_class XX, ZZ, w_temp1, w_temp2, w, s_temp, s, ss, R, RR, XaddR, B_temp1, XXaddRR, B, ww, B2, h, Ylh, Y_temp1, Yrh;
	Point_proj result;
	sqr_fp(p.X, mod, &XX);
	sqr_fp(p.Z, mod, &ZZ);

	mul_fp(a, ZZ, mod, &w_temp1);
	mul_fp(3, XX, mod, &w_temp2);
	add_fp(w_temp1, w_temp2, mod, &w);

	mul_fp(2, p.Y, mod, &s_temp);
	mul_fp(s_temp, p.Z, mod, &s);

	sqr_fp(s, mod, &ss);

	mul_fp(s, ss, mod, &result.Z);

	mul_fp(p.Y, s, mod, &R);

	sqr_fp(R, mod, &RR);

	add_fp(p.X, R, mod, &XaddR);
	sqr_fp(XaddR, mod, &B_temp1);
	add_fp(XX, RR, mod, &XXaddRR);
	sub_fp(B_temp1, XXaddRR, mod, &B);

	sqr_fp(w, mod, &ww);
	mul_fp(2, B, mod, &B2);
	sub_fp(ww, B2, mod, &h);

	mul_fp(h, s, mod, &result.X);

	sub_fp(B, h, mod, &Y_temp1);
	mul_fp(w, Y_temp1, mod, &Ylh);
	mul_fp(2, RR, mod, &Yrh);
	sub_fp(Ylh, Yrh, mod, &result.Y);
	
	result.inf = false;
	return result;
}

//点の加算
//参考(https://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-add-1998-cmo-2)
Point_proj proj_add(const Point_proj& p, const Point_proj& q, const mpz_class& a, const mpz_class& mod)
{
	if (p.X == 0 && p.Z == 0) return q;
	if (q.X == 0 && q.Z == 0) return p;
	mpz_class t0, t1, u0, u1, t, tt, u, u2, v, w, u3, X3, Y3, Z3, w_temp1, w_temp2, addu;
	mpz_class Ytemp1, Ytemp2, Ytemp3, Ytemp4;
	Point_proj result;
	mul_fp(p.Y, q.Z, mod, &t0);
	mul_fp(q.Y, p.Z, mod, &t1);
	mul_fp(p.X, q.Z, mod, &u0);
	mul_fp(q.X, p.Z, mod, &u1);
	if(u0==u1 && t0==t1) return ec_double_proj(p, a, mod);	//2倍算
	sub_fp(t0, t1, mod, &t);
	sub_fp(u0, u1, mod, &u);
	mul_fp(u, u, mod, &u2);
	mul_fp(p.Z, q.Z, mod, &v);
	sqr_fp(t, mod, &tt);
	mul_fp(tt, v, mod, &w_temp1);
	add_fp(u0, u1, mod, &addu);
	mul_fp(u2, addu, mod, &w_temp2);
	sub_fp(w_temp1, w_temp2, mod, &w);
	mul_fp(u, u2, mod, &u3);
	
	mul_fp(u, w, mod, &result.X);

	mul_fp(u0, u2, mod, &Ytemp1);
	sub_fp(Ytemp1, w, mod, &Ytemp2);
	mul_fp(Ytemp2, t, mod, &Ytemp3);
	mul_fp(t0, u3, mod, &Ytemp4);
	sub_fp(Ytemp3, Ytemp4, mod, &result.Y);

	mul_fp(u3, v, mod, &result.Z);
	result.inf = false;
	return result;
}
