use std::f64::consts::PI;

use iir_design as lib;

use lib::elliptics::jacobi;

use lib::domain::{FilterShapeKind, PassBandFilterKind};

fn main() {
    todo!()
}

//		private void Design(int order, double delta, double scale, double eps)
fn design(
    order: u32,
    delta: f64,
    scale: f64,
    eps: f64,
    pass_band_kind: PassBandFilterKind,
    filter_shape: FilterShapeKind,
    f_sample: f64, // sampling frequency
    fc_low: f64,   // cut off frequency (low side)
    fc_high: f64,  // cut off frequency (high side)
) {
    let mut v0: f64 = 0.;
    let mut kk = 0.;
    let mut k = 0.;
    let mut k2 = 0.;

    //			double fa, ang, wc, cgam, cbp;

    //if (passBand_.IsChecked("LPF", "HPF")) {
    //    ang = fcLow_ * PI / fSample_;
    //    fa = fSample_ / 2.0;
    //} else {
    //    ang = Abs(fcHigh_ - fcLow_) * PI / fSample_;
    //    fa = fcHigh_;
    //}
    let (angle, fa) = match pass_band_kind {
        PassBandFilterKind::LowPass | PassBandFilterKind::HighPass => {
            (fc_low * PI / f_sample, f_sample / 2.0)
        }
        PassBandFilterKind::BandPass | PassBandFilterKind::BandRejection => {
            ((fc_high - fc_low).abs() * PI / f_sample, fc_high)
        }
    };

    let (cgam, kk, k, v0, k2, wc, cbp) = match filter_shape {
        FilterShapeKind::Elliptic => {
            //if (!EllipInp(order, fa, eps, ang, out cgam, out Kk, out k, out v0,
            //  out k2, out wc, out cbp)) return;
            todo!()
        }
        _ => {
            // 周波数変換（LPF → BPF）
            let wc = angle.tan();

            let cgam = (PI * (fa + fc_low) / f_sample).cos() / (angle.cos());
            let fa = 2.0 * PI * fc_high / f_sample;
            let cbp = (cgam - fa.cos()) / fa.sin();
            (cgam, kk, k, v0, k2, wc, cbp)
        }
    };

    //sPlane(order, delta, k2, Kk, v0, k, wc, out np, out nz);
    let (np, nz) = s_plane(order, delta, k2, kk, v0, k, wc);
    //			AddZero(Bilinear(np, nz, cgam, Tan(ang), cbp, wc));
    //			ToDirect(cgam, scale);
    //			ToCascade();
}

//		// λ平面における位置を捜す
//		private void LambdaPlane(double order, double eps, double wr,
//								 out double Kk, out double v0, out double k2, out double wc)
fn lambda_plane(order: f64, eps: f64, wr: f64) -> (f64, f64, f64, f64) {
    let wc = 1.0;
    let k: f64 = wc / wr;
    let k2 = k * k;
    let kk = jacobi::ellip_c(1.0 - k2);
    let k1d: f64 = jacobi::modulus(order * jacobi::ellip_c(k2) / kk); // k1'
    let v0 = kk * jacobi::ellip_i(1.0 / eps, k1d) / (order * jacobi::ellip_c(1.0 - k1d * k1d)); // 式(7.82)
    (kk, v0, k2, wc)
}

//		// sButterworth, sChebyshev で使用する関数
//		private int sButChe(int order, int np, double a, double b)
//		{
//			int nz = 0;
//			double theta = ((order & 1) == 1) ? 0.0 : PI/(2.0*order);
//
//			for (int n=0; n<np; n++)	// 極
//			{
//				s_pz_[n] = new Complex(-a*Cos(theta), b*Sin(theta));
//				if (formType_.IsChecked("逆")) s_pz_[n] = 1.0/s_pz_[n];		// s → 1/s
//				theta += PI/order;
//			}
//
//			if (formType_.IsChecked("逆"))
//			{
//				for (int n=0; n<(order+1)/2; n++)
//					s_pz_[n+np] = ImaginaryOne/Cos((2*n+1)*PI/(2.0*order));
//				nz = order/2;
//			}
//
//			if (passBand_.IsChecked("HPF", "BRF"))	// HPF, BRF の場合
//			{
//				for (int n=0; n<np; n++) s_pz_[n] = 1.0/s_pz_[n]; // s → 1/s
//				nz = np;
//				if (passBand_.IsChecked("BRF")) nz += order/2;
//
//				if (formType_.IsChecked("バ", "チ"))
//					for (int j=0; j<nz; j++) s_pz_[j+np] = 0.0;
//				else
//					for (int j=0; j<nz; j++) s_pz_[j+np] = 1.0/s_pz_[j+np];
//			}
//
//			return nz;
//		}

//		// Butterworth フィルタ：s 平面の極と零点の計算（遮断角周波数 = 1）
//		private int sButterworth(int order, int np)
//		{
//			return (sButChe(order, np, 1.0, 1.0));
//		}

//		// Chebyshev フィルタ：s 平面の極と零点の計算（遮断角周波数 = 1）
//		private int sChebyshev(int order, double delta, int np)
//		{
//			double x = Sqrt(delta*delta - 1.0);					// 逆Chebyshev
//			if (formType_.IsChecked("チ")) x = 1.0/x;			// Chebyshev
//			double alpha = Log(x + Sqrt(x*x  + 1.0))/order;		// arcsinh(x)
//			return sButChe(order, np, Sinh(alpha), Cosh(alpha));
//		}

//		// Elliptic フィルタ：s 平面の極と零点の計算（遮断角周波数 = 1）
//		private int sElliptic(int order, double m, double Kk, double v0, double k, double wc, int np)
//		{
//			double sn1, cn1, dn1, sn, cn, dn;
//
//			int nz = order/2;
//			Jacobi.SnCnDn(v0, m, out sn1, out cn1, out dn1);	// 式(7.85) sn', cn', dn' を求める
//			for (int n=0; n<nz; n++)		// 零点
//			{
//				Jacobi.SnCnDn((order-1-2*n)*Kk/order, 1.0-m, out sn, out cn, out dn);
//				s_pz_[n+np] = new Complex(0.0, wc/(k*sn));
//			}
//			for (int n=0; n<np; n++)		// 極
//			{
//				Jacobi.SnCnDn((order-1-2*n)*Kk/order, 1.0-m, out sn, out cn, out dn);
//				double r = k*sn*sn1;
//				double denom = cn1*cn1 + r*r;
//				s_pz_[n] = new Complex(-wc*cn*dn*sn1*cn1/denom, wc*sn*dn1/denom);
//			}
//			if (passBand_.IsChecked("HPF", "BRF"))
//			{
//				for (int n=0; n<(np+nz); n++) s_pz_[n] = 1.0/s_pz_[n];
//				while (np > nz) s_pz_[np+nz++] = 0.0;
//			}
//			return (nz);
//		}

// s 平面の極と零点の計算（遮断角周波数 = 1）
//		private void sPlane(int order, double delta, double m, double Kk, double v0, double k, double wc,
//							out int np, out int nz)
fn s_plane(order: u32, delta: f64, m: f64, kk: f64, v0: f64, k: f64, wc: f64) -> (u32, u32) {
    //			s_pz_ = new Complex[(zOrd_+1)*2];
    //			Array.Clear(s_pz_, 0, (zOrd_+1)*2);
    let np = (order + 1) / 2;
    let nz = 0;
    //			switch (formType_.Checked)
    //			{
    //				case "バ": nz = sButterworth(order, np); break;
    //				case "チ": nz = sChebyshev(order, delta, np); break;
    //				case "逆": nz = sChebyshev(order, delta, np); break;
    //				case "連": nz = sElliptic(order, m, Kk, v0, k, wc, np); break;
    //			}
    (np, nz)
}

//		// 双一次ｚ変換
//		private int Bilinear(int np, int nz, double cgam, double c, double cbp, double wc)
//		{
//			const double ZERO = 1.0E-16;
//
//			z_pz_ = new Complex[zOrd_*2];
//			Array.Clear(z_pz_, 0, z_pz_.Length);
//
//			double xc = formType_.IsChecked("連") ? c : wc;
//			int nc = np;
//			int jt = -1;
//			int j = 0;
//
//			int nLoop = (nz == 0) ? 1 : 2;
//			for (int i_pz=0; i_pz<nLoop; i_pz++)	// s 平面から z 平面への変換
//			{
//				for (int n=0; n<nc; n++)
//				{
//					Complex r = s_pz_[j++];
//					if (passBand_.IsChecked("LPF", "HPF"))
//					{								// LPF, HPF の場合
//						z_pz_[++jt] = (1.0 + xc*r)/(1.0 - xc*r);
//						if (Abs(r.Imaginary) >= ZERO)
//						{
//							jt++;
//							z_pz_[jt] = Conjugate(z_pz_[jt-1]);
//						}
//					}
//					else
//					{								// BPF, BRF の場合
//						Complex cx = r*(formType_.IsChecked("チ") ? cbp : c);
//						Complex ca = 1.0 - cx;
//						Complex cb = -2.0*cgam;
//						Complex b4ac = Sqrt(cb*cb - 4.0*(1.0 - cx*cx));
//						Complex root = (b4ac - cb)/(2.0*ca);
//
//						if (jt >= zOrd_*2-1) break;
//						z_pz_[++jt] = root;
//
//						if (Abs(root.Imaginary) > ZERO) z_pz_[++jt] = Conjugate(root);
//
//						if ((Abs(r.Imaginary) > ZERO) || (Abs(root.Imaginary) < ZERO))
//						{
//							root = -(b4ac + cb)/(2.0*ca);
//							z_pz_[++jt] = root;
//							if (Abs(root.Imaginary) > ZERO) z_pz_[++jt] = Conjugate(root);
//						}
//					}
//				}
//				nc = nz;
//			}
//			return jt;
//		}

//		// 残りの零点を付加
//		private void AddZero(int jt)
pub fn add_zero(
    jt: usize,
    filter_shape: FilterShapeKind,
    z_ord: usize,
    pass_band: PassBandFilterKind,
) {
    //			if (!formType_.IsChecked("連"))		// Elliptic 以外
    match filter_shape {
        FilterShapeKind::Elliptic => {
            while 2 * z_ord - 1 > jt {
                //					z_pz_[++jt] = -1.0;
                //if (passBand_.IsChecked("BPF", "BRF")) {
                //  z_pz_[++jt] = 1.0;
                //}
                todo!()
            }
        }
        _ => {
            //				while (2*zOrd_ - 1 > jt)
            //				{
            //					if (!passBand_.IsChecked("HPF")) z_pz_[++jt] = -1.0;
            //					if (passBand_.IsChecked("HPF", "BPF")) z_pz_[++jt] = 1.0;
            //				}
            todo!()
        }
    }
    //			else								// Elliptic
    todo!()
}

//		// 直接形の係数に変換
//		private void ToDirect(double cgam, double scale)
pub fn to_direct(cgam: f64, scale: f64, z_ord: usize) {
    //			double an, bn;
    //
    //let ak_ = new double[zOrd_+1];
    let ak = vec![0.; z_ord + 1];
    let bk = vec![0.; z_ord + 1];
    //			Expand(ak_, 0);
    //			Expand(bk_, zOrd_);
    let a: f64 = 1.0;

    //			if (passBand_.IsChecked("BPF"))	// BPF の場合
    //			{
    //				double gam = PI/2.0 - Asin(cgam);
    //				int mh = zOrd_/2;
    //				an = ak_[mh];
    //				bn = bk_[mh];
    //				double ai = 0.0;
    //				if (mh > ((zOrd_/4)*2))
    //				{
    //					ai = 1.0;
    //					an = 0.0;
    //					bn = 0.0;
    //				}
    //				for (int j=1; j<=mh; j++)
    //				{
    //					a = gam*j - ai*PI/2.0;
    //					int jh = mh + j;
    //					int jl = mh - j;
    //					an += Cos(a)*(ak_[jh] + (1.0 - 2.0*ai)*ak_[jl]);
    //					bn += Cos(a)*(bk_[jh] + (1.0 - 2.0*ai)*bk_[jl]);
    //				}
    //			}
    //			else							// BPF 以外
    //			{
    //				if (passBand_.IsChecked("HPF")) a = -1.0;
    //				an = 1.0;
    //				bn = 1.0;
    //				for (int j=1; j<=zOrd_; j++)
    //				{
    //					an = a*an + ak_[j];
    //					bn = a*bn + bk_[j];
    //				}
    //			}
    //
    //			double gain = an/(bn*scale);
    //			if (!formType_.IsChecked("連") && (bn == 0)) gain = 1.0;
    //
    //			for (int j=0; j<=zOrd_; j++) bk_[j] = gain*bk_[j];
    todo!()
}

//		// 多項式の係数に展開
//		private void Expand(double[] ck, int order)
//		{
//			Complex[] xk = new Complex[zOrd_+1];
//
//			xk[0] = 1.0;
//			for (int n=1; n<=zOrd_; n++) xk[n] = 0.0;
//			for (int n=0; n<zOrd_; n++)
//				for (int k=0; k<=n; k++)
//					xk[n-k+1] = xk[n-k+1] - z_pz_[n+order]*xk[n-k];
//			for (int n=0; n<=zOrd_; n++) ck[n] = xk[n].Real;
//		}

//		// 縦続形の係数に変換
//		private void ToCascade()
//		{
//			cm_ = new Coefs[(zOrd_+1)/2];
//
//			if (formType_.IsChecked("逆") && ((zOrd_ % 2) == 1) &&
//				(passBand_.IsChecked("LPF") || passBand_.IsChecked("HPF")))
//			{
//				Complex tmp = z_pz_[2*zOrd_-1];
//				for (int n=2*zOrd_-1; n>zOrd_; n--) z_pz_[n] = z_pz_[n-1];
//				z_pz_[zOrd_] = tmp;
//			}
//
//			int m = 0;
//			for (int j=0; j<zOrd_; j++)
//			{
//				int jz = j + zOrd_;
//				if (Abs(z_pz_[j].Imaginary) > 1.0E-16)				// 複素極のステージ
//				{
//					cm_[m].a1 = 2.0*z_pz_[j].Real;					// a1m
//					double abs_a2m = Abs(z_pz_[j]);
//					cm_[m].a2 = -abs_a2m*abs_a2m; 					// a2m
//					cm_[m].b0 = 1.0;								// b0m
//					cm_[m].b1 = -(z_pz_[jz] + z_pz_[jz+1]).Real;	// b1m
//					cm_[m].b2 = (z_pz_[jz]*z_pz_[jz+1]).Real;		// b2m
//					j++;
//				}
//				else												// 実極のステージ
//				{
//					if ((zOrd_ % 2) != 0)							// 奇数次の場合
//					{
//						cm_[m].a1 = z_pz_[j].Real;					// a1m
//						cm_[m].a2 = 0.0;							// a2m = 0
//						cm_[m].b0 = 1.0;							// b0m
//						cm_[m].b1 = -z_pz_[jz].Real;				// b1m
//						cm_[m].b2 = 0.0;							// b2m = 0
//					}
//
//					else											// 偶数次の場合
//					{
//						cm_[m].a1 = (z_pz_[j] + z_pz_[j+1]).Real;	// a1m
//						cm_[m].a2 = -(z_pz_[j]*z_pz_[j+1]).Real;	// a2m
//						cm_[m].b0 = 1.0;							// b0m
//						cm_[m].b1 = -(z_pz_[jz] + z_pz_[jz+1]).Real;// b1m
//						cm_[m].b2 = (z_pz_[jz]*z_pz_[jz+1]).Real;	// b2m
//						j++;
//					}
//				}
//				m++;
//			}
//
//			gain_ = Gain();		// 利得定数を縦続形の係数から求める
//
//			// ソート，極の絶対値が小さい順，ただし実極は先頭になる
//			cm_ = cm_.OrderBy(x => Abs(x.a2)).ToArray();
//		}

//		// 利得定数を縦続形の係数から求める
//		private double Gain()
//		{
//			double g = 1.0;
//			switch (formType_.Checked)
//			{
//				case "バ": g = Sqrt(2.0); break;
//				case "チ": g = Pow(10, dbRrpl_/20.0); break;
//				case "逆": g = Pow(10, dbD_/20.0); break;
//				case "連": g = Pow(10, dbRrpl_/20.0); break;
//			}
//			return 1.0/AbsHw(fcLow_*2.0*PI/fSample_, g);
//		}

//		// L∞ノルムによるスケーリング
//		private void Scaling()
//		{
//			double[] norm = new double[(zOrd_+1)/2];
//			scaleC_ = new double[(zOrd_+1)/2+1];
//			cmS_ = new Coefs[(zOrd_+1)/2];
//
//			for (int n=0; n<(zOrd_+1)/2; n++) norm[n] = PeakPicking(n);
//
//			scaleC_[0] = 1.0/norm[0];
//			for (int n=1; n<(zOrd_+1)/2; n++) scaleC_[n] = norm[n-1]/norm[n];
//			scaleC_[(zOrd_+1)/2] = gain_*norm[(zOrd_+1)/2-1];
//
//			for (int n=0; n<cm_.Length; n++)
//			{
//				cmS_[n].a1 = cm_[n].a1;
//				cmS_[n].a2 = cm_[n].a2;
//				cmS_[n].b0 = scaleC_[n+1]*cm_[n].b0;
//				cmS_[n].b1 = scaleC_[n+1]*cm_[n].b1;
//				cmS_[n].b2 = scaleC_[n+1]*cm_[n].b2;
//			}
//		}

//		// 振幅特性の最大値検出
//		private double PeakPicking(int sct)
//		{
//			double amp_max = 0.0;
//			double dw = PI/1.0E4;
//
//			for (double w=0; w<=PI; w+=dw)
//				amp_max = Max(amp_max, Hw_n(w, sct));
//			return amp_max;
//		}

//		// 入力点からオーバーフローに注目する点までの伝達関数
//		private double Hw_n(double w, int sect)
//		{
//			Complex z = Exp(-w*ImaginaryOne);	// 1/z
//			Complex z2 = z*z;
//
//			Complex h = 1.0/(1.0 - cm_[0].a1*z - cm_[0].a2*z2); // 第１セクション
//			if (sect == 0) return Abs(h);
//
//			for (int n=1; n<=sect; n++) 						// 第２セクション以降
//				h = h*(cm_[n-1].b0 + cm_[n-1].b1*z + cm_[n-1].b2*z2)
//					/(1.0 - cm_[n].a1*z - cm_[n].a2*z2);
//			return Abs(h);
//		}
//	}

// Elliptic フィルタの設計パラメータ入力その他
//private bool EllipInp(double order, double fa, double eps, double arg,
//					  out double cgam, out double Kk, out double k, out double u,
//					  out double k2, out double wc, out double cbp)
fn ellip_inp(
    order: f64,
    fa: f64,
    eps: f64,
    arg: f64,
    fc_low: f64,
    f_sample: f64,
    db_d: f64,
    pass_band_kind: PassBandFilterKind,
)
//					  out double cgam, out double Kk, out double k, out double u,
//					  out double k2, out double wc, out double cbp)
{
    //	wc = u = Kk = k2 = k = cbp = 0;

    let _cgam = ((fa + fc_low) * PI / f_sample).cos() / arg.cos();
    //	if (!dBdownInput()) return false;

    let m1: f64 = eps / (10.0f64.powf(db_d / 10.0) - 1.0).sqrt();
    let m1 = m1 * m1;
    let k = jacobi::modulus(jacobi::ellip_c(m1) / (order as f64 * jacobi::ellip_c(1.0 - m1)));

    let wr: f64 = match pass_band_kind {
        PassBandFilterKind::HighPass | PassBandFilterKind::BandRejection => k,
        PassBandFilterKind::LowPass | PassBandFilterKind::BandPass => 1. / k,
    };

    let f3: f64 = match pass_band_kind {
        PassBandFilterKind::LowPass | PassBandFilterKind::HighPass => {
            //		f3 = Atan(Tan(arg)*wr)*fSample_/PI;		// LPF, HPF の場合
            todo!()
        }
        PassBandFilterKind::BandPass | PassBandFilterKind::BandRejection => {
            //	{											// BPF, BRF の場合
            //		double wrtan = wr*Tan(arg);
            //		double wrtan2 = wrtan*wrtan;
            //		double tmp = (cgam + Sqrt(wrtan2*((1.0 - cgam*cgam) + wrtan2)))/(1.0 + wrtan2);
            //		f3 = (PI/2.0 - Asin(tmp))*fSample_/(2.0*PI);
            //	}
            todo!()
        }
    };

    let angle = f3 * PI / f_sample;
    let cang = angle.cos();
    let sang = angle.sin();

    //	if (passBand_.IsChecked("LPF", "HPF"))
    //		wr = sang/(cang*Tan(arg));						// LPF, HPF の場合
    //	else												// BPF, BRF の場合
    //	{
    //		double q = cang*cang - sang*sang;
    //		sang = 2.0*cang*sang;
    //		cang = q;
    //		wr = (cgam - cang)/(sang*Tan(arg));
    //	}
    let wr: f64 = match pass_band_kind {
        PassBandFilterKind::LowPass | PassBandFilterKind::HighPass => todo!(),
        PassBandFilterKind::BandPass | PassBandFilterKind::BandRejection => todo!(),
    };

    //	if (passBand_.IsChecked("HPF", "BRF")) {
    //wr = 1.0/wr;	// HPF, BRF の場合
    //}
    let wr = wr.abs();
    let cbp = wr;
    //LambdaPlane(order, eps, wr, out Kk, out u, out k2, out wc);
    let (kk, v0, k2, wc) = lambda_plane(order, eps, wr);

    //	return true;
}

// 阻止域の減衰量の入力
//private bool dBdownInput()
fn db_down_input() -> f64 {
    //bool
    //			is_plus = textBox4.ToDouble(out dbD_, v => v > 0.0,
    //									 "阻止域の減衰量は正の値を入力してください．");
    //return is_plus
    todo!()
}
