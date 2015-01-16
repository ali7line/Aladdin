
typedef union
#ifdef __cplusplus
	YYSTYPE
#endif
 {
	SYMBOL     *sym;
	Inst   	  *inst;
	int	   narg;
} YYSTYPE;
extern YYSTYPE yylval;
# define QUIT 257
# define END_OF_FILE 258
# define BREAK 259
# define SET_UNITS_OFF 260
# define SET_UNITS_ON 261
# define SECTION 262
# define MATERIAL 263
# define TYPE 264
# define FIBER 265
# define WHILE 266
# define IF 267
# define ELSE 268
# define THEN 269
# define FOR 270
# define NUMBER 271
# define STRING 272
# define PRINT 273
# define UNDEF 274
# define BLTINVar_MATRIX 275
# define BLTINVar_QUANTITY 276
# define BLTIN_MATRIX 277
# define BLTIN1_MATRIX 278
# define BLTIN2_MATRIX 279
# define BLTIN3_MATRIX 280
# define BLTIN_QUANTITY 281
# define BLTIN1_QUANTITY 282
# define BLTIN2_QUANTITY 283
# define BLTIN3_QUANTITY 284
# define NODE_QUANT 285
# define MESH 286
# define LINK_NODE 287
# define ADD_ELMT 288
# define SECT_ATTR 289
# define ELMT_ATTR 290
# define MATL_ATTR 291
# define FIB_ATTR 292
# define UNIT 293
# define MAP 294
# define LDOF 295
# define TO 296
# define GDOF 297
# define BLTIN_FE_FUNC 298
# define BLTIN1_FE_FUNC 299
# define VAR 300
# define QUAN 301
# define VECT 302
# define MATX 303
# define DIMENSION 304
# define LENGTH 305
# define MASS 306
# define GTIME 307
# define OR 308
# define AND 309
# define GT 310
# define GE 311
# define LT 312
# define LE 313
# define EQ 314
# define NE 315
# define UNARYPLUS 316
# define UNARYMINUS 317
# define NOT 318
