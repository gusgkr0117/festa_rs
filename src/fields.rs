pub mod FpFESTA {
    const N: usize = 21;
    const BITLEN: usize = 1293;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x43FFFFFFFFFFFFFF,
        0x9AC3245C4D7BE6B3,
        0x21D7DCD797059B7B,
        0x8F19A73E323F6569,
        0x841FED4773CFDB16,
        0x02979D50DD13D09A,
        0x01712922BAF59934,
        0xBD1C756E54F72C15,
        0xF6B3CF47C54370FE,
        0xCEC87BD4C1480F2B,
        0x11CF13E54B11406F,
        0x000000000000176C,
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0xA200000000000000,
        0xCD61922E26BDF359,
        0x90EBEE6BCB82CDBD,
        0x478CD39F191FB2B4,
        0x420FF6A3B9E7ED8B,
        0x014BCEA86E89E84D,
        0x80B894915D7ACC9A,
        0x5E8E3AB72A7B960A,
        0xFB59E7A3E2A1B87F,
        0xE7643DEA60A40795,
        0x08E789F2A588A037,
        0x0000000000000BB6,
    ];
    const R_VAL: [u64; N] = [
        0x000AEE091BECF8C7,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x2400000000000000,
        0x78FB641CC3C350C6,
        0xA7023CB1847A0623,
        0x81CA1B3E72AA42FF,
        0x375C61B4D97A8C2B,
        0x36365530C09FB122,
        0xC85A4827D5B4CC8C,
        0x07F34CE9CA1014FC,
        0x5515312398D806A2,
        0x1051CF6B62E881A8,
        0xD49418553ECD1FE6,
        0x0000000000000A37,
    ];
    const MINUS_R_VAL: [u64; N] = [
        0xFFF511F6E4130738,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x1FFFFFFFFFFFFFFF,
        0x21C7C03F89B895ED,
        0x7AD5A026128B9558,
        0x0D4F8BFFBF952269,
        0x4CC38B929A554EEB,
        0xCC6148201C741F78,
        0x3916E0FAE540CCA7,
        0xB52928848AE71718,
        0xA19E9E242C6B6A5C,
        0xBE76AC695E5F8D83,
        0x3D3AFB900C442089,
        0x0000000000000D34,
    ];
    const DR_VAL: [u64; N] = [
        0x0015DC1237D9F18E,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x4800000000000000,
        0xF1F6C8398786A18C,
        0x4E04796308F40C46,
        0x0394367CE55485FF,
        0x6EB8C369B2F51857,
        0x6C6CAA61813F6244,
        0x90B4904FAB699918,
        0x0FE699D3942029F9,
        0xAA2A624731B00D44,
        0x20A39ED6C5D10350,
        0xA92830AA7D9A3FCC,
        0x000000000000146F,
    ];
    const TR_VAL: [u64; N] = [
        0x0020CA1B53C6EA56,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x2800000000000000,
        0xD02F07F9FDCE0B9F,
        0xD32ED93CF66876EE,
        0xF644AA7D25BF6395,
        0x21F537D7189FC96B,
        0xA00B624164CB42CC,
        0x579DAF54C628CC70,
        0x5ABD714F093912E1,
        0x088BC4230544A2E7,
        0x622CF26D677175CD,
        0x6BED351A71561F42,
        0x000000000000073B,
    ];
    const QR_VAL: [u64; N] = [
        0x002BB8246FB3E31D,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x4C00000000000000,
        0x492A6C16C1915C65,
        0x7A3115EE7AE27D12,
        0x780EC5BB9869A695,
        0x5951998BF21A5597,
        0xD641B772256AF3EE,
        0x1FF7F77C9BDD98FC,
        0x62B0BE38D34927DE,
        0x5DA0F5469E1CA989,
        0x727EC1D8CA59F775,
        0x40814D6FB0233F28,
        0x0000000000001173,
    ];
    const R2_VAL: [u64; N] = [
        0x9F42DE4F788045C2,
        0x1C2EB545BBE1E6CA,
        0x2C7698FD862CBB36,
        0x00FDEE42F690BB02,
        0x37B26D78444A90CA,
        0x95622499BE0F470F,
        0xC62D0EAF8B7DC9CC,
        0xF242919A892AA2DF,
        0x601AD6354D69F3CB,
        0x074F84AEB6FA3F2B,
        0xF14BEC8FE0ABEB2D,
        0xD8339C56ACAD0170,
        0xF06FF606EE809357,
        0x88D41B5189622D47,
        0xB6C213020828FB73,
        0x9F8125EC514BBAAC,
        0xAC872B8F8CFA8AD7,
        0x711147A40E5E32CA,
        0xBCF167C3A3CBB9FE,
        0x76F7DA793EAA61EF,
        0x0000000000000D90,
    ];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [
        0x57CFDEF044149D6F,
        0xCE41D70241B51263,
        0xFB982E57F54C7E28,
        0xA4679409D1E483C7,
        0x011DC2D5A0BAADE8,
        0x8AAA62A138C96A66,
        0xD222A8115270BA21,
        0x06C5848819FD1F77,
        0x9BB2C052D2153FD6,
        0x066E9B66267C36CA,
        0x3D5097316A500C8B,
        0xC61ADE5C138AE9CA,
        0x3C7F105387EA0D08,
        0x9179D77EB3F9D1DA,
        0x94B6937BF7B9E519,
        0xDF141B79C2C7EFA0,
        0xE30A188AF2ECCFE1,
        0xF7E4F400401E6695,
        0x7C3B040888C99764,
        0xDCF405E68822CC43,
        0x00000000000005D2,
    ];
    const TDEC_VAL: [u64; N] = [
        0x1C2EB545BBE1E6CA,
        0x2C7698FD862CBB36,
        0x00FDEE42F690BB02,
        0x37B26D78444A90CA,
        0x95622499BE0F470F,
        0xC62D0EAF8B7DC9CC,
        0xF242919A892AA2DF,
        0x601AD6354D69F3CB,
        0x8F4F84AEB6FA3F2B,
        0x71325BC65B41105A,
        0x79727C59FFDA3258,
        0xFF3C908B8ED06518,
        0xE1D183D9987FBC43,
        0x43072C6DA68FE17C,
        0x8A0380573C4B2FAD,
        0x108DB7C73303ED54,
        0xD97B27C377C7B8EA,
        0x792E18D4CD60A8A6,
        0x628D06B8251703F4,
        0x3D43319001A8F373,
        0x0000000000000E92,
    ];
    const WIN_LEN: usize = 5;
    const SQRT_EH: [u8; 133] = [
        4, 26, 12, 13, 30, 29, 21, 9, 28, 2, 9, 6, 12, 13, 14, 15, 27, 12, 1, 14, 25, 11, 19, 27,
        23, 14, 8, 18, 22, 18, 29, 7, 18, 17, 15, 14, 26, 12, 28, 17, 22, 24, 22, 31, 28, 25, 29,
        8, 13, 31, 7, 8, 8, 13, 2, 26, 19, 8, 23, 1, 21, 14, 30, 18, 2, 0, 13, 18, 25, 26, 11, 23,
        2, 9, 10, 2, 23, 0, 20, 2, 12, 25, 29, 9, 5, 23, 21, 14, 28, 8, 15, 29, 15, 24, 13, 8, 5,
        30, 17, 30, 28, 25, 26, 30, 11, 25, 3, 16, 20, 0, 19, 26, 27, 3, 18, 29, 28, 23, 1, 8, 17,
        24, 18, 10, 30, 9, 28, 25, 17, 0, 27, 14, 1,
    ];
    const SQRT_EL: usize = 126;
    const FOURTH_ROOT_EH: [u8; 132] = [
        2, 13, 22, 6, 31, 30, 26, 4, 14, 17, 4, 3, 22, 6, 23, 23, 13, 22, 0, 23, 28, 21, 25, 29,
        11, 7, 4, 9, 11, 25, 30, 3, 25, 24, 7, 7, 13, 6, 30, 8, 11, 12, 27, 15, 30, 28, 14, 20, 22,
        31, 3, 4, 20, 6, 1, 29, 9, 20, 27, 16, 10, 7, 15, 9, 1, 16, 6, 25, 12, 29, 21, 11, 17, 4,
        5, 17, 11, 0, 10, 1, 22, 28, 30, 20, 18, 27, 10, 7, 14, 20, 23, 30, 7, 28, 6, 20, 2, 31, 8,
        15, 30, 12, 13, 31, 21, 28, 1, 8, 10, 16, 9, 29, 29, 1, 25, 14, 30, 27, 0, 20, 8, 12, 9, 5,
        31, 4, 30, 28, 8, 16, 13, 23,
    ];
    const FOURTH_ROOT_EL: usize = 126;
    const P1: u64 = 3143667320;
    const P1DIV_M: u64 = 6755719943463976019;

    crate::finitefield::fp_gen::define_fp_core! {}

    #[cfg(test)]
    mod tests {
        crate::finitefield::fp_gen::define_fp_tests! {}
    }
}

pub mod FpFESTAExt {
    use num_bigint::BigUint;

    use super::FpFESTA::Fp;
    const NQR_RE: Fp = Fp::new([
        0x145E4F94B005484D,
        0x437E234B4E5CF5D6,
        0x5D29F56B6150B469,
        0xD1DCEFD038BFD33D,
        0xC3C772943D654936,
        0xCEAC95C36C49B427,
        0xE282B9FB7D9CBE13,
        0x0677E9713A992BF4,
        0xB7ECEF2567197DCC,
        0x8E9451C2845DFB09,
        0xA1229BB35705E1B9,
        0x38BC621169F0A3C7,
        0x215B414B1AB53288,
        0x1C3384A030CC972B,
        0xD7E9BB1787954CE7,
        0xCD570AD9A9444C6F,
        0x3EC28136E459A065,
        0x96D89B4337340AC7,
        0xA92B8E9DEBBBEB2B,
        0xA4D7F84010F736DF,
        0x0000000000001047,
    ]);

    /// Torsion basis of the elliptic curve E0 : y^2 = x^3 + 6 * x^2 + x
    pub const P0x: Fp = Fp::new([
        0xB2DAFEECE8954A26,
        0x7C30CFF9C9295BE9,
        0x9BA18C569D14412C,
        0xE0E3B63074BEC0AE,
        0x398044D58E03B86D,
        0xC8D0AAC8160AC9FF,
        0xC91CC769AB73DA65,
        0x7A4D84690DD347CC,
        0x0E86C1850765979C,
        0xD04F4C29371AB193,
        0x57567A513DE96E94,
        0x201F334E7B86DB61,
        0xFA905410682435F3,
        0xCF94E4C18F3CC5A3,
        0xF871351240711D39,
        0x1259D6CBBC4FF5F9,
        0xFA608AA257D9E807,
        0x0298E79660009EED,
        0x1EADA2117D2996A5,
        0xCF14C3787CD20B3E,
        0x0000000000000B22,
    ]);
    pub const P0y: Fp = Fp::new([
        0xA926641D45B8A118,
        0xDA4B2092E262BD3B,
        0x33169DAA5116A1BE,
        0xE5FBD874E2E731E2,
        0x30E54C3ED0E8784D,
        0x26DF07E9EDB68209,
        0xCCCFB6B687006A44,
        0x1EADD94C491708A5,
        0x9A1D9825932BEC3F,
        0x7C22FD699B440B45,
        0x690DE842938D6D48,
        0xAF6EFDB934F17EC4,
        0x181A55AA14F763E2,
        0x23D08570DDCAE0B5,
        0x905F8561A4F79ACA,
        0xFD89660855A08BDA,
        0x2E5D73B8E567A03C,
        0x8046C1FDD942F3F6,
        0x090140DD8B71CD57,
        0xF3C01B38A431D9B5,
        0x0000000000000522,
    ]);
    pub const Q0x: Fp = Fp::new([
        0x54FCDF4057795A20,
        0x9CE149FD474BB801,
        0xCA9F07D54A0CD609,
        0xEC42C1D4B79FE4E0,
        0xF9EC70EC0214422E,
        0xF58BECE0A9404F35,
        0xB46EE18B995B528C,
        0x09516EA0D30F9713,
        0xCDCAA700A0ECA732,
        0x31CCD53267606A86,
        0x653751C72E982B44,
        0xB3455F058CD7E1E1,
        0x503C9CCE4755C2FC,
        0x194E42014A0E15F2,
        0xA54BA2F12BE6825A,
        0xCC5CF225A610EB25,
        0x60022E3792FA436B,
        0xDFDA10165CC7A209,
        0x42567069E2CA09ED,
        0x8FF030C55A5C6DD6,
        0x0000000000001123,
    ]);
    pub const Q0y: Fp = Fp::new([
        0x1CAE4655B7317ED4,
        0x2CDA0F3338024BF0,
        0x4009114682F7EFFE,
        0xABF26A97D0B7307D,
        0x5FC9EFAA795D975B,
        0x069D2E0D2D3C5BF0,
        0x16239137BE42D449,
        0xC6E44780FA7EE4A7,
        0x1995B5152D368056,
        0x1E3D6B1970DB5F9D,
        0xB34DC5908C4ABE6E,
        0x8FA82803240603FE,
        0x4A746A5A1FD93A43,
        0xAE5D9F859CD30A11,
        0xA49D4B643D8B6D87,
        0xA611BAD0D873B395,
        0x1EDD64CAE0736797,
        0xE0D4298F428D2ED7,
        0x8717FB8AEE0C371E,
        0xDC3012D1E48DD904,
        0x000000000000138F,
    ]);

    /// Public parameters
    pub const D1: [u32; 9] = [
        0xECF246B9, 0x6A652D42, 0x20A0E9BD, 0x8F6346C8, 0x62F334BB, 0x9769376C, 0x3CA861E7,
        0x6C1B23A8, 0x00000001,
    ];
    pub const D1_FACTORED: [(u32, u32); 16] = [
        (3, 6),
        (19, 2),
        (29, 2),
        (37, 2),
        (83, 2),
        (139, 2),
        (167, 2),
        (251, 2),
        (419, 2),
        (421, 2),
        (701, 2),
        (839, 2),
        (1009, 2),
        (1259, 2),
        (3061, 2),
        (3779, 2),
    ];
    pub const D2: [u32; 9] = [
        0xE9035DDF, 0x9EDC26E4, 0xB3FDD8F0, 0x0B756A33, 0x69CB451C, 0x6FA41A34, 0x4D255546,
        0xB6B52466, 0x0000000D,
    ];
    pub const D2_FACTORED: [(u32, u32); 18] = [
        (5, 4),
        (7, 3),
        (11, 2),
        (13, 2),
        (17, 2),
        (41, 2),
        (43, 2),
        (71, 2),
        (89, 2),
        (127, 2),
        (211, 2),
        (281, 2),
        (503, 2),
        (631, 2),
        (2309, 2),
        (2521, 2),
        (2647, 2),
        (2729, 2),
    ];

    pub const DA: [u32; 9] = [
        0x76F4A401, 0x5FB7F347, 0x9A503BB9, 0x2EC47996, 0x6E83E49F, 0xE7CC1C9D, 0xEC0B3EF0,
        0xA85D5DC4, 0x00012BF3,
    ];
    pub const DA_FACTORED: [(u32, u32); 11] = [
        (59, 2),
        (3023, 2),
        (3359, 2),
        (4409, 2),
        (5039, 2),
        (6299, 2),
        (6719, 2),
        (9181, 2),
        (19531, 2),
        (22679, 2),
        (41161, 2),
    ];

    pub const DA1: [u32; 3] = [0x92CE3979, 0xCB60C53C, 0x01B2BF4B];
    pub const DA1_FACTORED: [(u32, u32); 4] = [(59, 1), (6299, 1), (6719, 1), (9181, 1)];
    pub const DA2: [u32; 6] = [
        0x71F8A4C9, 0x5E8DBA36, 0x32DE49A5, 0x7AEF9314, 0x016200DD, 0x00B0A040,
    ];
    pub const DA2_FACTORED: [(u32, u32); 7] = [
        (3023, 1),
        (3359, 1),
        (4409, 1),
        (5039, 1),
        (19531, 1),
        (22679, 1),
        (41161, 1),
    ];

    pub const M1: [u32; 5] = [0x6FC3B771, 0xCA5A2369, 0xB11FABF6, 0x232ED4F9, 0x00001121];
    pub const M2: [u32; 3] = [0x3A8A70CD, 0x69735CD2, 0x52C628EA];

    /// original 632 plus extra 2bits torsion information
    /// for the fast theta based isogeny operation
    pub const L_POWER: u32 = 634;
    /// length of the theta based (2,2)^l-isogeny
    pub const THETA_ISOGENY_LENGTH: u32 = 632;

    /// order of the torsion basis
    pub const BASIS_ORDER: [u32; 41] = [
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x44000000, 0x4D7BE6B3,
        0x9AC3245C, 0x97059B7B, 0x21D7DCD7, 0x323F6569, 0x8F19A73E, 0x73CFDB16, 0x841FED47,
        0xDD13D09A, 0x02979D50, 0xBAF59934, 0x01712922, 0x54F72C15, 0xBD1C756E, 0xC54370FE,
        0xF6B3CF47, 0xC1480F2B, 0xCEC87BD4, 0x4B11406F, 0x11CF13E5, 0x0000176C,
    ];

    pub const BIGUINT_MODULUS: [u32; 41] = [
        0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
        0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
        0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x43FFFFFF, 0x4D7BE6B3,
        0x9AC3245C, 0x97059B7B, 0x21D7DCD7, 0x323F6569, 0x8F19A73E, 0x73CFDB16, 0x841FED47,
        0xDD13D09A, 0x02979D50, 0xBAF59934, 0x01712922, 0x54F72C15, 0xBD1C756E, 0xC54370FE,
        0xF6B3CF47, 0xC1480F2B, 0xCEC87BD4, 0x4B11406F, 0x11CF13E5, 0x0000176C,
    ];

    /// Precomputed (2,2)-isogeny strategy
    pub const THETA_STRATEGY: [usize; 631] = [
        233, 144, 89, 55, 42, 34, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1,
        1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1,
        1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 13, 8, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8,
        5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1,
        1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3,
        2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 55, 34, 21, 13, 8, 5, 3, 2, 1,
        1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3,
        2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1,
        1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1,
        1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 89, 55, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2,
        1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1,
        1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1,
        1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1,
        3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1,
        1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1,
        1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2,
        1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1,
        1, 1, 1, 1, 2, 1, 1, 1,
    ];

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        // crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}
