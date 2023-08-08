pub enum PassBandFilterKind {
    LowPass,
    HighPass,
    BandPass,
    BandRejection,
}

pub enum FilterShapeKind {
    Butterworth,
    Chebyshev,
    //sChebyshev
    Elliptic,
}
