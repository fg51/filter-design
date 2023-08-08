//	Jacobi の楕円関数に関連するメソッド

pub mod jacobi {

    use std::f64::consts::PI;

    //	public static class Jacobi
    //	{
    //		// K/K' から母数（modulus）を求める
    //		public static double Modulus(double u)
    pub fn modulus(u: f64) -> f64 {
        //			double q = Exp(-PI*u);
        let q: f64 = (-PI * u).exp();
        let mut a: f64 = 1.0;
        let mut b: f64 = 1.0;
        let mut c: f64 = 1.0;
        let mut d: f64 = q;

        let count: usize = 15;
        for _ in 0..count {
            a += 2.0 * c * d;
            c *= d * d;
            b += c;
            d *= q;
            if c < 1.0E-16 {
                break;
            }
        }
        //    if c >= 1.0E-16 {
        //MessageBox.Show("Modulus: 収束しません．");
        //    }
        let a = b / a;
        4.0 * (q.sqrt()) * a * a
    }

    // 第１種完全楕円積分
    // 　算術幾何平均による方法を利用
    //	 山内他編：“電子計算機のための数値計算法Ⅲ”，p.258，培風館，1972 参照
    //	 1-x*x に対する値を計算
    //		public static double EllipC(double x)
    pub fn ellip_c(x: f64) -> f64 {
        let mut a = 1.0;
        let mut b = x.sqrt();

        let count: usize = 20;
        for _n in 0..count {
            let at = (a + b) / 2.0;
            b = (a * b).sqrt();
            a = at;
            if (a - b) / a < 1.0E-8 {
                break;
            }
        }
        //			if (n >= COUNT)
        //			{
        //				MessageBox.Show("EllipC: 収束しません．");
        //				return 0;
        //			}

        (PI / 2.0) / a
    }

    // 第１種不完全楕円積分
    //	 ただし，kc = sqrt(1-k*k)
    //		public static double EllipI(double x, double kc)
    pub fn ellip_i(x: f64, kc: f64) -> f64 {
        let mut a: f64 = 1.0;
        let mut b: f64 = kc;
        let mut y: f64 = 1.0 / x;
        let mut i: f64 = 0.;

        //			int n;
        let count: usize = 15;
        for _ in 0..count {
            let bt = a * b;
            a += b;
            b = 2.0 * bt.sqrt();
            y -= bt / y;
            if y == 0.0 {
                y = bt.sqrt() * 1.0E-10;
            }
            if (a - b).abs() < a * 1.0E-10 {
                break;
            }
            i *= 2.;
            if y < 0.0 {
                i += 1.;
            }
        }
        //			if (n >= COUNT)
        //			{
        //				MessageBox.Show("EllipI: 収束しません．");
        //				return 0;
        //			}
        if y < 0.0 {
            i += 1.
        }
        ((a / y).atan() + PI * i) / a
    }

    // Jacobi の楕円関数，sn(u,kc), cn(u,kc), dn(u,kc) の計算
    //	 ただし，kc = 1-k*k
    //		public static bool SnCnDn(double u, double kc, out double sn, out double cn, out double dn)
    pub fn sn_cn_dn(u: f64, kc: f64) //, out double sn, out double cn, out double dn)
    {
        //			sn = 0.0;
        //			cn = 1.0;
        //			dn = 1.0;
        //			if (u == 0.0) return true;

        let mut a = 1.0;
        let mut b = kc.sqrt();
        const COUNT: usize = 16;
        let mut aa = [0.; COUNT];
        let mut bb = [0.; COUNT];
        //			int n;
        let mut n = 0;
        for i in 0..COUNT {
            aa[i] = a;
            bb[i] = b;
            let at = (a + b) / 2.0;
            b = (a * b).sqrt();
            a = at;
            n = i;
            if (a - b) / a < 1.0E-8 {
                break;
            }
        }
        //			if (n >= COUNT)
        //			{
        //				MessageBox.Show("SnCnDn: 収束しません．");
        //				return false;
        //			}

        let mut c: f64 = a / (u * a).tan();
        let mut d = 1.0;
        let mut k = n as i32;
        while k >= 0 {
            //			for (int k=n; k>=0; k--)
            let e = c * c / a;
            c *= d;
            a = aa[k as usize];
            d = (e + bb[k as usize]) / (e + a);
            k -= 1;
        }
        let sn = 1.0 / (1.0 + c * c).sqrt();
        let _cn = sn * c;
        let _dn = d;
        //			return true;
        //		}
        todo!()
    }
}
