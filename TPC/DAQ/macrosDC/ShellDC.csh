
#filename = $1
root.exe -l << EOF
.L Fun4All_TPC_UnpackDC.C
//Fun4All_TPC_UnpackPRDF(200000,"$1");
//Fun4All_TPC_UnpackDC(20000,"$1");
//Fun4All_TPC_UnpackDC(4000,"$1");
Fun4All_TPC_UnpackDC(3000000,"$1");
//Fun4All_TPC_UnpackDC(450,"$1");
.q
EOF

