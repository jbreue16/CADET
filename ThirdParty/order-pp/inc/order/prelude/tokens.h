// (C) Copyright Vesa Karvonen 2004.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE.)

#define ORDER_PP_DEF_8has_tuples ORDER_PP_FN_CM(1,8HAS_TUPLES,0IS_ANY)
#define ORDER_PP_8HAS_TUPLES(P,t,...) (,ORDER_PP_IS_EDIBLE(,P##t)(,8true,8false),P##__VA_ARGS__)

#define ORDER_PP_DEF_8drop_tuples ORDER_PP_FN_CM(1,8DROP_TUPLES,0IS_ANY)
#define ORDER_PP_8DROP_TUPLES(P,ts,...) (,ORDER_PP_TUPLES_TERMINATE(ORDER_PP_DROP_TUPLES_A ts##P),P##__VA_ARGS__)
#define ORDER_PP_DROP_TUPLES_A(...) ORDER_PP_DROP_TUPLES_B
#define ORDER_PP_DROP_TUPLES_B(...) ORDER_PP_DROP_TUPLES_A
#define ORDER_PP_ORDER_PP_DROP_TUPLES_A
#define ORDER_PP_ORDER_PP_DROP_TUPLES_B

#define ORDER_PP_DEF_8take_tuples ORDER_PP_FN_CM(1,8TAKE_TUPLES,0IS_ANY)
#define ORDER_PP_8TAKE_TUPLES(P,ts,...) (,ORDER_PP_SCAN(ORDER_PP_EAT ORDER_PP_LPAREN ORDER_PP_TAKE_TUPLES_A ts##P ORDER_PP_RPAREN),P##__VA_ARGS__)
#define ORDER_PP_TAKE_TUPLES_A(...) )((__VA_ARGS__))ORDER_PP_EAT ORDER_PP_LPAREN ORDER_PP_TAKE_TUPLES_B
#define ORDER_PP_TAKE_TUPLES_B(...) )((__VA_ARGS__))ORDER_PP_EAT ORDER_PP_LPAREN ORDER_PP_TAKE_TUPLES_A

#define ORDER_PP_DEF_8tokens_to_seq ORDER_PP_FN_CM(1,8TOKENS_TO_SEQ,0IS_ANY)
#define ORDER_PP_8TOKENS_TO_SEQ(P,ts,...) (,ORDER_PP_FX(PSEQ_TO_SEQ,(,ORDER_PP_FX(PSEQ_REVERSE,(,ORDER_PP_PM(,ORDER_PP_IS_EDIBLE(,P##ts)(ORDER_PP_TOKENS_TO_SEQ_,EDIBLE,INEDIBLE)(,P##ts,ORDER_PP_TOKEN),ORDER_PP_TOKEN,))))),P##__VA_ARGS__)

#define ORDER_PP_DEF_8tokens_to_seq_with ORDER_PP_FN_CM(2,8TOKENS_TO_SEQ_WITH,0IS_ANY,0IS_ANY)
#define ORDER_PP_8TOKENS_TO_SEQ_WITH(P,ts,m,...) (,ORDER_PP_FX(PSEQ_TO_SEQ,(,ORDER_PP_FX(PSEQ_REVERSE,(,ORDER_PP_PM(,ORDER_PP_IS_EDIBLE(,P##ts)(ORDER_PP_TOKENS_TO_SEQ_,EDIBLE,INEDIBLE)(,P##ts,P##m),P##m,))))),P##__VA_ARGS__)

#define ORDER_PP_8TOKENS_TO_SEQ_LOOP(P,x,ts,stop_nil,m,...) (,ORDER_PP_IS_EDIBLE(,P##ts)(ORDER_PP_TOKENS_TO_SEQ_,EDIBLE,INEDIBLE)(,P##ts,P##m),P##m,)(,P##x)

#define ORDER_PP_TOKENS_TO_SEQ_EDIBLE(P,ts,m) ORDER_PP_TOKENS_TO_SEQ_EDIBLE_TAKE ts##P,
#define ORDER_PP_TOKENS_TO_SEQ_EDIBLE_TAKE(x) (x),8TOKENS_TO_SEQ_LOOP,

#define ORDER_PP_TOKENS_TO_SEQ_INEDIBLE(P,ts,m) ORDER_PP_FX(TOKENS_TO_SEQ_INEDIBLE_TAKE,m##_##ts),8STOP_NIL
#define ORDER_PP_TOKENS_TO_SEQ_INEDIBLE_TAKE(x) x,8TOKENS_TO_SEQ_LOOP,

#define ORDER_PP_TOKEN ORDER_PP_TOKEN

#define ORDER_PP_TOKEN_0 (0)
#define ORDER_PP_TOKEN_1 (1)
#define ORDER_PP_TOKEN_2 (2)
#define ORDER_PP_TOKEN_3 (3)
#define ORDER_PP_TOKEN_4 (4)
#define ORDER_PP_TOKEN_5 (5)
#define ORDER_PP_TOKEN_6 (6)
#define ORDER_PP_TOKEN_7 (7)
#define ORDER_PP_TOKEN_8 (8)
#define ORDER_PP_TOKEN_9 (9)

// Detail

#define ORDER_PP_TUPLES_TERMINATE(x) ORDER_PP_PRIMITIVE_CAT(ORDER_PP_,x)