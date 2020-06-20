{
    "index": "world", 
    "drawstyle": "wireframe", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "material": "air", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "", 
    "invisible": 1, 
    "size": [
        2000.0, 
        2000.0, 
        2000.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "dark_box", 
    "drawstyle": "wireframe", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "material": "cardboard", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "world", 
    "size": [
        1250.0, 
        1250.0, 
        1250.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "air_box", 
    "drawstyle": "wireframe", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "material": "air", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "dark_box", 
    "size": [
        1220.0, 
        1220.0, 
        1220.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "acrylic_wall", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "r_min": 497.5, 
    "material": "acrylic_uvt", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 522.9, 
    "size_z": 625.0, 
    "mother": "air_box", 
    "position": [
        0.0, 
        0.0, 
        150.0
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "cd", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "material": "water", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 497.5, 
    "size_z": 625.0, 
    "mother": "air_box", 
    "position": [
        0.0, 
        0.0, 
        150.0
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "teflon_wall", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.3
    ], 
    "r_min": 481.825, 
    "material": "teflon_black", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 485.0, 
    "size_z": 625.0, 
    "mother": "cd", 
    "position": [
        0.0, 
        0.0, 
        0.0
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "acrylic_bottom", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "material": "acrylic_uvt", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 522.9, 
    "size_z": 12.7, 
    "mother": "air_box", 
    "position": [
        0.0, 
        0.0, 
        -487.70000000000005
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "acrylic_flange", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "r_min": 522.9, 
    "material": "acrylic_uvt", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 575.0, 
    "size_z": 12.7, 
    "mother": "air_box", 
    "position": [
        0.0, 
        0.0, 
        762.3
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "acrylic_lid", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "material": "acrylic_uvt", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 575.0, 
    "size_z": 12.7, 
    "mother": "air_box", 
    "position": [
        0.0, 
        0.0, 
        787.7
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "center_hole", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "material": "water", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 88.9, 
    "size_z": 12.7, 
    "mother": "acrylic_lid", 
    "position": [
        0.0, 
        0.0, 
        0.0
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "side_hole", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "material": "water", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 88.9, 
    "size_z": 12.7, 
    "mother": "acrylic_lid", 
    "position": [
        -293.7, 
        0.0, 
        0.0
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "center_port_wall", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "r_min": 88.9, 
    "material": "acrylic_uvt", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 101.6, 
    "size_z": 60.0, 
    "mother": "air_box", 
    "position": [
        0.0, 
        0.0, 
        860.4
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "side_port_wall", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "r_min": 88.9, 
    "material": "acrylic_uvt", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 101.6, 
    "size_z": 60.0, 
    "mother": "air_box", 
    "position": [
        -293.7, 
        0.0, 
        860.4
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "center_port_flange", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "r_min": 101.6, 
    "material": "acrylic_uvt", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 138.1, 
    "size_z": 12.7, 
    "mother": "air_box", 
    "position": [
        0.0, 
        0.0, 
        907.6999999999999
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "side_port_flange", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "r_min": 101.6, 
    "material": "acrylic_uvt", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 138.1, 
    "size_z": 12.7, 
    "mother": "air_box", 
    "position": [
        -293.7, 
        0.0, 
        907.6999999999999
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "center_port_lid", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "material": "acrylic_uvt", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 138.1, 
    "size_z": 12.7, 
    "mother": "air_box", 
    "position": [
        0.0, 
        0.0, 
        933.1
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "side_port_lid", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "material": "acrylic_uvt", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 138.1, 
    "size_z": 12.7, 
    "mother": "air_box", 
    "position": [
        -293.7, 
        0.0, 
        933.1
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "center_overflow", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "material": "water", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 88.9, 
    "size_z": 30.0, 
    "mother": "air_box", 
    "position": [
        0.0, 
        0.0, 
        830.4
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "side_overflow", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        0.0, 
        0.5
    ], 
    "material": "water", 
    "valid_begin": [
        0, 
        0
    ], 
    "r_max": 88.9, 
    "size_z": 30.0, 
    "mother": "air_box", 
    "position": [
        -293.7, 
        0.0, 
        830.4
    ], 
    "type": "tube", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "optical_cookies_0",
    "name": "GEO",
    "color": [
        1.0,
        1.0,
        0.0,
        0.5
    ],
    "material": "opticalcookies",
    "valid_begin": [
        0,
        0
    ],
    "r_max": 25.4,
    "size_z": 1.0,
    "mother": "air_box",
    "position": [
        381.35,
        -30.79375,
        -501.40000000000005
    ],
    "type": "tube",
    "valid_end": [
        0,
        0
    ]
}
{
    "index": "optical_cookies_1",
    "name": "GEO",
    "color": [
        1.0,
        1.0,
        0.0,
        0.5
    ],
    "material": "opticalcookies",
    "valid_begin": [
        0,
        0
    ],
    "r_max": 25.4,
    "size_z": 1.0,
    "mother": "air_box",
    "position": [
        131.35,
        -30.79375,
        -501.40000000000005
    ],
    "type": "tube",
    "valid_end": [
        0,
        0
    ]
}
{
    "index": "optical_cookies_2",
    "name": "GEO",
    "color": [
        1.0,
        1.0,
        0.0,
        0.5
    ],
    "material": "opticalcookies",
    "valid_begin": [
        0,
        0
    ],
    "r_max": 25.4,
    "size_z": 1.0,
    "mother": "air_box",
    "position": [
        -118.65,
        -30.79375,
        -501.40000000000005
    ],
    "type": "tube",
    "valid_end": [
        0,
        0
    ]
}
{
    "index": "optical_cookies_3",
    "name": "GEO",
    "color": [
        1.0,
        1.0,
        0.0,
        0.5
    ],
    "material": "opticalcookies",
    "valid_begin": [
        0,
        0
    ],
    "r_max": 25.4,
    "size_z": 1.0,
    "mother": "air_box",
    "position": [
        -368.65,
        -30.79375,
        -501.40000000000005
    ],
    "type": "tube",
    "valid_end": [
        0,
        0
    ]
}
{
    "index": "optical_cookies_4",
    "name": "GEO",
    "color": [
        1.0,
        1.0,
        0.0,
        0.5
    ],
    "material": "opticalcookies",
    "valid_begin": [
        0,
        0
    ],
    "r_max": 25.4,
    "size_z": 1.0,
    "mother": "air_box",
    "position": [
        -165.025,
        280.375,
        801.4
    ],
    "type": "tube",
    "valid_end": [
        0,
        0
    ]
}
{
    "index": "optical_cookies_5",
    "name": "GEO",
    "color": [
        1.0,
        1.0,
        0.0,
        0.5
    ],
    "material": "opticalcookies",
    "valid_begin": [
        0,
        0
    ],
    "r_max": 25.4,
    "size_z": 1.0,
    "mother": "air_box",
    "position": [
        -165.025,
        -128.425,
        801.4
    ],
    "type": "tube",
    "valid_end": [
        0,
        0
    ]
}
{
    "index": "optical_cookies_6",
    "name": "GEO",
    "color": [
        1.0,
        1.0,
        0.0,
        0.5
    ],
    "material": "opticalcookies",
    "valid_begin": [
        0,
        0
    ],
    "r_max": 25.4,
    "size_z": 1.0,
    "mother": "air_box",
    "position": [
        381.35,
        169.20625,
        -501.40000000000005
    ],
    "type": "tube",
    "valid_end": [
        0,
        0
    ]
}
{
    "index": "optical_cookies_7",
    "name": "GEO",
    "color": [
        1.0,
        1.0,
        0.0,
        0.5
    ],
    "material": "opticalcookies",
    "valid_begin": [
        0,
        0
    ],
    "r_max": 25.4,
    "size_z": 1.0,
    "mother": "air_box",
    "position": [
        131.35,
        169.20625,
        -501.40000000000005
    ],
    "type": "tube",
    "valid_end": [
        0,
        0
    ]
}
{
    "index": "pmts_S0", 
    "pmt_detector_type": "idpmt", 
    "enable": 1, 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "sensitive_detector": "/mydet/pmt/inner", 
    "material": null, 
    "pmt_model": "r7723_S0", 
    "pos_table": "PMTINFO_S0", 
    "valid_begin": [
        0, 
        0
    ], 
    "orientation": "manual", // "point" means all PMTs at a point
	     	             // "manual" means there is a dir_x, dir_y, dir_z
	                     // in pos_table for manual orientation of PMTs
    "mother": "air_box", 
    "efficiency_correction": 1.0, 
    "type": "pmtarray", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "pmts_S1", 
    "pmt_detector_type": "idpmt", 
    "enable": 1, 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "sensitive_detector": "/mydet/pmt/inner", 
    "material": null, 
    "pmt_model": "r7723_S1", 
    "pos_table": "PMTINFO_S1", 
    "valid_begin": [
        0, 
        0
    ], 
    "orientation": "manual", // "point" means all PMTs at a point
	     	             // "manual" means there is a dir_x, dir_y, dir_z
	                     // in pos_table for manual orientation of PMTs
    "mother": "air_box", 
    "efficiency_correction": 1.0, 
    "type": "pmtarray", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "pmts_S2", 
    "pmt_detector_type": "idpmt", 
    "enable": 1, 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "sensitive_detector": "/mydet/pmt/inner", 
    "material": null, 
    "pmt_model": "r7723_S2", 
    "pos_table": "PMTINFO_S2", 
    "valid_begin": [
        0, 
        0
    ], 
    "orientation": "manual", // "point" means all PMTs at a point
	     	             // "manual" means there is a dir_x, dir_y, dir_z
	                     // in pos_table for manual orientation of PMTs
    "mother": "air_box", 
    "efficiency_correction": 1.0, 
    "type": "pmtarray", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "pmts_S3", 
    "pmt_detector_type": "idpmt", 
    "enable": 1, 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "sensitive_detector": "/mydet/pmt/inner", 
    "material": null, 
    "pmt_model": "r7723_S3", 
    "pos_table": "PMTINFO_S3", 
    "valid_begin": [
        0, 
        0
    ], 
    "orientation": "manual", // "point" means all PMTs at a point
	     	             // "manual" means there is a dir_x, dir_y, dir_z
	                     // in pos_table for manual orientation of PMTs
    "mother": "air_box", 
    "efficiency_correction": 1.0, 
    "type": "pmtarray", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "pmts_S4", 
    "pmt_detector_type": "idpmt", 
    "enable": 1, 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "sensitive_detector": "/mydet/pmt/inner", 
    "material": null, 
    "pmt_model": "r7723_S4", 
    "pos_table": "PMTINFO_S4", 
    "valid_begin": [
        0, 
        0
    ], 
    "orientation": "manual", // "point" means all PMTs at a point
	     	             // "manual" means there is a dir_x, dir_y, dir_z
	                     // in pos_table for manual orientation of PMTs
    "mother": "air_box", 
    "efficiency_correction": 1.0, 
    "type": "pmtarray", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "pmts_S5", 
    "pmt_detector_type": "idpmt", 
    "enable": 1, 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "sensitive_detector": "/mydet/pmt/inner", 
    "material": null, 
    "pmt_model": "r7723_S5", 
    "pos_table": "PMTINFO_S5", 
    "valid_begin": [
        0, 
        0
    ], 
    "orientation": "manual", // "point" means all PMTs at a point
	     	             // "manual" means there is a dir_x, dir_y, dir_z
	                     // in pos_table for manual orientation of PMTs
    "mother": "air_box", 
    "efficiency_correction": 1.0, 
    "type": "pmtarray", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "pmts_S6", 
    "pmt_detector_type": "idpmt", 
    "enable": 1, 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "sensitive_detector": "/mydet/pmt/inner", 
    "material": null, 
    "pmt_model": "r7723_S6", 
    "pos_table": "PMTINFO_S6", 
    "valid_begin": [
        0, 
        0
    ], 
    "orientation": "manual", // "point" means all PMTs at a point
	     	             // "manual" means there is a dir_x, dir_y, dir_z
	                     // in pos_table for manual orientation of PMTs
    "mother": "air_box", 
    "efficiency_correction": 1.0, 
    "type": "pmtarray", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "pmts_S7", 
    "pmt_detector_type": "idpmt", 
    "enable": 1, 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "sensitive_detector": "/mydet/pmt/inner", 
    "material": null, 
    "pmt_model": "r7723_S7", 
    "pos_table": "PMTINFO_S7", 
    "valid_begin": [
        0, 
        0
    ], 
    "orientation": "manual", // "point" means all PMTs at a point
	     	             // "manual" means there is a dir_x, dir_y, dir_z
	                     // in pos_table for manual orientation of PMTs
    "mother": "air_box", 
    "efficiency_correction": 1.0, 
    "type": "pmtarray", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_platform", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "air_box", 
    "position": [
        0.0, 
        0.0, 
        -509.9
    ], 
    "size": [
        650.0, 
        650.0, 
        9.5
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_bar_group", 
    "name": "GEO", 
    "color": [
        0.0, 
        0.0, 
        0.0, 
        0.2
    ], 
    "material": "air", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_platform", 
    "size": [
        631.0, 
        631.0, 
        9.5
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_bottom", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "air_box", 
    "position": [
        0.0, 
        0.0, 
        -1190.9
    ], 
    "size": [
        509.5, 
        509.5, 
        9.5
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_bottom_cut", 
    "name": "GEO", 
    "color": [
        0.0, 
        0.0, 
        0.0, 
        0.2
    ], 
    "material": "air", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_bottom", 
    "size": [
        490.5, 
        490.5, 
        9.5
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_bottom_bar", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_bottom_cut", 
    "size": [
        9.5, 
        500.0, 
        9.5
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_left_group", 
    "drawstyle": "wireframe", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "material": "air", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "air_box", 
    "invisible": 1, 
    "position": [
        -500.0, 
        0.0, 
        -850.4
    ], 
    "size": [
        9.5, 
        650.0, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_center_group", 
    "drawstyle": "wireframe", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "material": "air", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "air_box", 
    "invisible": 1, 
    "position": [
        0.0, 
        0.0, 
        -850.4
    ], 
    "size": [
        9.5, 
        650.0, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_right_group", 
    "drawstyle": "wireframe", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        0.5
    ], 
    "material": "air", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "air_box", 
    "invisible": 1, 
    "position": [
        500.0, 
        0.0, 
        -850.4
    ], 
    "size": [
        9.5, 
        650.0, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_bar_0", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_bar_group", 
    "position": [
        0.0, 
        -500.0, 
        0.0
    ], 
    "size": [
        631.0, 
        9.5, 
        9.5
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_bar_1", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_bar_group", 
    "position": [
        0.0, 
        -300.0, 
        0.0
    ], 
    "size": [
        631.0, 
        9.5, 
        9.5
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_bar_2", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_bar_group", 
    "position": [
        0.0, 
        -100.0, 
        0.0
    ], 
    "size": [
        631.0, 
        9.5, 
        9.5
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_bar_3", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_bar_group", 
    "position": [
        0.0, 
        100.0, 
        0.0
    ], 
    "size": [
        631.0, 
        9.5, 
        9.5
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_bar_4", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_bar_group", 
    "position": [
        0.0, 
        300.0, 
        0.0
    ], 
    "size": [
        631.0, 
        9.5, 
        9.5
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_bar_5", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_bar_group", 
    "position": [
        0.0, 
        500.0, 
        0.0
    ], 
    "size": [
        631.0, 
        9.5, 
        9.5
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_left_0", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_left_group", 
    "position": [
        0.0, 
        -500.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_left_1", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_left_group", 
    "position": [
        0.0, 
        -300.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_left_2", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_left_group", 
    "position": [
        0.0, 
        -100.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_left_3", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_left_group", 
    "position": [
        0.0, 
        100.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_left_4", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_left_group", 
    "position": [
        0.0, 
        300.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_left_5", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_left_group", 
    "position": [
        0.0, 
        500.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_center_0", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_center_group", 
    "position": [
        0.0, 
        -500.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_center_1", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_center_group", 
    "position": [
        0.0, 
        -300.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_center_2", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_center_group", 
    "position": [
        0.0, 
        -100.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_center_3", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_center_group", 
    "position": [
        0.0, 
        100.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_center_4", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_center_group", 
    "position": [
        0.0, 
        300.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_center_5", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_center_group", 
    "position": [
        0.0, 
        500.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_right_0", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_right_group", 
    "position": [
        0.0, 
        -500.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_right_1", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_right_group", 
    "position": [
        0.0, 
        -300.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_right_2", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_right_group", 
    "position": [
        0.0, 
        -100.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_right_3", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_right_group", 
    "position": [
        0.0, 
        100.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_right_4", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_right_group", 
    "position": [
        0.0, 
        300.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "frame_leg_right_5", 
    "name": "GEO", 
    "color": [
        0.0, 
        1.0, 
        1.0, 
        0.2
    ], 
    "material": "aluminum", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "frame_leg_right_group", 
    "position": [
        0.0, 
        500.0, 
        0.0
    ], 
    "size": [
        9.5, 
        9.5, 
        331.0
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "scintillation_strip_0", 
    "sensitive_detector": "/mydet/veto/genericchamber", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        1.0
    ], 
    "material": "pvc", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "air_box", 
    "position": [
        266.70000000000005, 
        24.41575, 
        1133.5375
    ], 
    "size": [
        50.8, 
        44.45, 
        3.175
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "scintillation_strip_1", 
    "sensitive_detector": "/mydet/veto/genericchamber", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        1.0
    ], 
    "material": "pvc", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "air_box", 
    "position": [
        267.04999999999995, 
        24.47925, 
        887.4376
    ], 
    "size": [
        50.8, 
        57.15, 
        4.7244
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "scintillation_strip_2", 
    "sensitive_detector": "/mydet/veto/genericchamber", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        1.0
    ], 
    "material": "pvc", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "air_box", 
    "position": [
        319.08750000000003, 
        26.193749999999998, 
        1144.65
    ], 
    "size": [
        50.8, 
        50.8, 
        3.175
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "scintillation_strip_3", 
    "sensitive_detector": "/mydet/veto/genericchamber", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        1.0
    ], 
    "material": "pvc", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "air_box", 
    "position": [
        273.05, 
        94.26575, 
        896.5435
    ], 
    "size": [
        50.8, 
        57.15, 
        4.3815
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "scintillation_strip_4", 
    "sensitive_detector": "/mydet/veto/genericchamber", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        1.0
    ], 
    "material": "pvc", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "air_box", 
    "position": [
        -174.625, 
        -141.7125, 
        -1123.8875
    ], 
    "size": [
        153.9875, 
        358.40625, 
        5.08
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
{
    "index": "scintillation_strip_5", 
    "sensitive_detector": "/mydet/veto/genericchamber", 
    "name": "GEO", 
    "color": [
        1.0, 
        1.0, 
        1.0, 
        1.0
    ], 
    "material": "pvc", 
    "valid_begin": [
        0, 
        0
    ], 
    "mother": "air_box", 
    "position": [
        179.3875, 
        -131.7625, 
        -1123.8875
    ], 
    "size": [
        152.4, 
        357.98125, 
        5.08
    ], 
    "type": "box", 
    "valid_end": [
        0, 
        0
    ]
}
