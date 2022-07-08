#include "fp.h"

//Weierstrass Curve
Point ec_add(const Point& p, const Point& q, const mpz_class& a, const mpz_class& mod)
{
	//cout << "p:" << p.x << " " << p.y << endl;
	//cout << "q:" << q.x << " " << q.y << endl;
	if (p.inf == 1) return q;
	if (q.inf == 1) return p;

	mpz_class lambda_temp, y_1, lambda,px2,temp;
	mpz_class lh, rh;
	Point result;


	//cout << "Px,Qx:" << p.x << "," << q.x << endl;
	if (p.x == q.x) {
		//P + (-P) = 0
		//cout << "Py,-Py:" << p.y << "," << -q.y << endl;
		if ((p.y + q.y == mod) || (p.y == q.y && p.y == 0)) {
			result.x = 0; result.y = 0;
			result.inf = true;
			return result;
		}
		//lambda=(3*x1*x1+a)/(2*y1)
		mul_fp(p.x, p.x, mod, &px2);
		mul_fp(3, px2, mod, &temp);
		add_fp(temp, a, mod, &lambda_temp);

		mul_fp(2, p.y, mod, &y_1);

		div_fp(lambda_temp, y_1, mod, &lambda);
	}
	else {
		//x1!=x2のとき，
		// cout << "Add" << endl;
		//lambda=(y2-y1)/(x2-x1)
		sub_fp(q.y, p.y, mod, &lh);
		sub_fp(q.x, p.x, mod, &rh);
		div_fp(lh, rh, mod, &lambda);
	}

	//P+Qのx座標
	mpz_class lambda2, x_temp;
	mul_fp(lambda, lambda, mod, &lambda2);
	add_fp(p.x, q.x, mod, &x_temp);
	sub_fp(lambda2, x_temp, mod, &result.x);

	//P+Qのy座標
	sub_fp(p.x, result.x, mod, &x_temp);

	mul_fp(lambda, x_temp, mod, &lambda2);
	sub_fp(lambda2, p.y, mod, &result.y);
	result.inf = false;
	return result;
}

Point copy_point(const Point& p) {
	Point result;
	result.x = p.x;
	result.y = p.y;
	result.inf = p.inf;
	return result;
}

bool sqrt_mod(const mpz_class& n, const mpz_class& p, mpz_class& result)
{
	mpz_class check;
	pow_fp(n, (p + 1) / 4, p, &result);

	pow_fp(result, 2, p, &check);
	if (check == n) {
		return true;
	}
	else return false;
}


/*** xを選んでからyを求める
	y^2=x^3+a*x+b ***/
Point gen_point_sqrt(const mpz_class& a, const mpz_class& b, const mpz_class& mod) {
	Point result;
	result.inf = false;
	mpz_class x, y, x_cube, x_sqr, y_temp, ax, rh, sqrt_result, testpow;
	result.x = random_fp(mod);
	//cout << "randX:" << result.x << endl;
	bool sqrt_check = 0;
	bool issqrt;

	while (sqrt_check == 0) {
		pow_fp(result.x, 3, mod, &x_cube);

		mul_fp(a, result.x, mod, &ax);

		add_fp(x_cube, ax, mod, &y_temp);
		add_fp(y_temp, b, mod, &rh);
		
		issqrt = sqrt_mod(rh, mod, sqrt_result);
		//sqrt_check = 1;
		if (issqrt == false) {	//整数じゃないとき
			result.x++;
			//cout << "rand:" << result.x << endl;
		}
		else {
			sqrt_check = 1;
			//mpz_sqrt(result.y.get_mpz_t(), rh.get_mpz_t());
			result.y = sqrt_result;
			//test
			pow_fp(result.y, 2, mod, &testpow);
			if (testpow != rh) cout << "違う" << endl;
			//cout <<"x,y :"<< result.x<<", " << result.y << endl;
		}
	}
	return result;
}

/*** y^2=x^3+a*x+bの点であれば0を返す ***/
size_t check_point(const Point& p, const mpz_class a, const mpz_class b, const mpz_class& mod)
{
	mpz_class x_cube, rh, lh, ax;
	pow_fp(p.x, 3, mod, &x_cube);

	mul_fp(a, p.x, mod, &ax);

	add_fp(x_cube, ax, mod, &rh);
	add_fp(rh, b, mod, &rh);

	pow_fp(p.y, 2, mod, &lh);

	if (lh == rh)	return 0;
	else			return 1;
}