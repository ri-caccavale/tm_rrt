%----------------------------------------------------------------------
% BEGIN LICENSE BLOCK
% Version: CMPL 1.1
%
% The contents of this file are subject to the Cisco-style Mozilla Public
% License Version 1.1 (the "License"); you may not use this file except
% in compliance with the License.  You may obtain a copy of the License
% at www.eclipse-clp.org/license.
% 
% Software distributed under the License is distributed on an "AS IS"
% basis, WITHOUT WARRANTY OF ANY KIND, either express or implied.  See
% the License for the specific language governing rights and limitations
% under the License. 
% 
% The Original Code is  The Gecode interface for ECLiPSe
% The Initial Developer of the Original Code is  Kish Shen
% Portions created by the Initial Developer are
% Copyright (C) 2009 Cisco Systems, Inc.  All Rights Reserved.
% 
% Contributor(s): Kish Shen
% 
% END LICENSE BLOCK
%----------------------------------------------------------------------
:- module(gfd).

:- lib(lists).
:- lib(constraint_pools).
:- lib(module_options).


:- import set_bip_error/1, get_bip_error/1 from sepia_kernel.

:- import
	timestamp_init/2,
        timestamp_update/2,
	timestamp_age/3,
        timestamp_older/4,
	request_fail_event/3 
   from sepia_kernel.

:- export op(700, xfx, [#::]).
:- export op(750, fx, [neg]).
:- export op(760, yfx, [and]).
:- export op(770, yfx, [or,xor]).
:- export op(790, yfx, [<=>]).
:- export op(780, yfx, [=>]).

:- export (#::)/2, (::)/2.
:- export (#::)/3, (::)/3.
:- export (#\=)/2, (#=)/2, (#<)/2, (#>)/2, (#>=)/2, (#=<)/2.
:- export (#\=)/3, (#=)/3, (#<)/3, (#>)/3, (#>=)/3, (#=<)/3.

:- export sumlist/2, sum/2, sum/3, sum/4, min/2, max/2, 
          scalar_product/4, scalar_product/5, divmod/4.
:- export all_le/2, all_lt/2, all_ge/2, all_gt/2, all_eq/2, all_ne/2.
:- export bool_channeling/3, inverse/2, inverse/4, inverse_g/2, inverse_g/4.
:- export ordered/2, lex_le/2, lex_lt/2, lex_ge/2, lex_gt/2, 
          lex_eq/2, lex_ne/2, mem/2, mem/3.
:- export alldifferent/1, alldifferent_cst/2, nvalues/3, gcc/2, occurrences/3, 
          atmost/3, atleast/3, count/4, count_matches/4, among/4. 
:- export element/3, element_g/3, precede/2, precede/3, 
          sorted/2, sorted/3, sorted_g/3. 
:- export circuit/1, circuit/3, circuit/4, circuit_offset/2, circuit_offset/4, 
          circuit_offset/5, circuit_g/1, circuit_g/3, circuit_g/4,
          circuit_offset_g/2,circuit_offset_g/4, circuit_offset_g/5,
          ham_path/3, ham_path/5, ham_path/6, ham_path_offset/4, 
          ham_path_offset/6, ham_path_offset/7, ham_path_g/3, ham_path_g/5, 
          ham_path_g/6, ham_path_offset_g/4, ham_path_offset_g/6, 
          ham_path_offset_g/7.
:- export disjunctive/2, disjunctive_optional/3,
          disjoint2/1, disjoint2_optional/1, 
          cumulative/4, cumulative_optional/5,
          cumulatives/5, cumulatives_min/5, 
          cumulatives_g/5, cumulatives_min_g/5. 
:- export sequence/5, sequence/4, bin_packing/3, bin_packing_g/3, bin_packing/4.
:- export table/2, table/3, extensional/4, regular/2.

:- export labeling/1, indomain/1, select_var/5, try_value/2.
%:- export labeling/3, indomain/2.
:- export is_in_domain/2, is_in_domain/3.
:- export search/6.
:- export gfd_update/0.

:- export (and)/2, (or)/2, (xor)/2, (<=>)/2, (=>)/2, neg/1.
:- export (and)/3, (or)/3, (xor)/3, (<=>)/3, (=>)/3, neg/2.

:- export get_min/2, get_max/2, get_median/2.
:- export get_bounds/3, get_integer_bounds/3, get_finite_integer_bounds/3.
:- export get_domain/2, get_domain_as_list/2, get_domain_size/2, 
          get_delta/2.
:- export get_constraints_number/2, get_weighted_degree/2.
:- export get_regret_lwb/2, get_regret_upb/2.
:- export impose_min/2, impose_max/2, impose_bounds/3, impose_domain/2,
          exclude/2, exclude_range/3.
:- export gfd_vars_impose_min/2, gfd_vars_impose_max/2, 
          gfd_vars_impose_bounds/3, gfd_vars_impose_domain/2,
          gfd_vars_exclude/2, gfd_vars_exclude_range/3, 
          gfd_vars_exclude_domain/2.
:- export is_solver_type/1, is_solver_var/1, is_exact_solver_var/1, integers/1.
:- export msg/3.

:- export gfd_maxint/1, gfd_minint/1.

:- export gfd_get_default/2, gfd_set_default/2.

:- tool('#\\='/2, '#\\=_body'/3).
:- tool('#='/2, '#=_body'/3).
:- tool('#<'/2, '#<_body'/3).
:- tool('#>'/2, '#>_body'/3).
:- tool('#>='/2, '#>=_body'/3).
:- tool('#=<'/2, '#=<_body'/3).

:- tool('#\\='/3, '#\\=_reif_body'/4).
:- tool('#='/3, '#=_reif_body'/4).
:- tool('#<'/3, '#<_reif_body'/4).
:- tool('#>'/3, '#>_reif_body'/4).
:- tool('#>='/3, '#>=_reif_body'/4).
:- tool('#=<'/3, '#=<_reif_body'/4).

:- tool('#\\=_c'/3, '#\\=_c'/4).
:- tool('#=_c'/3, '#=_c'/4).
:- tool('#<_c'/3, '#<_c'/4).
:- tool('#>_c'/3, '#>_c'/4).
:- tool('#>=_c'/3, '#>=_c'/4).
:- tool('#=<_c'/3, '#=<_c'/4).

:- tool('#\\=_reif_c'/4, '#\\=_reif_c'/5).
:- tool('#=_reif_c'/4, '#=_reif_c'/5).
:- tool('#<_reif_c'/4, '#<_reif_c'/5).
:- tool('#>_reif_c'/4, '#>_reif_c'/5).
:- tool('#>=_reif_c'/4, '#>=_reif_c'/5).
:- tool('#=<_reif_c'/4, '#=<_reif_c'/5).

/*
:- tool((=\=)/2, '#\\=_body'/3).
:- tool((=:=)/2, '#=_body'/3).
:- tool((<)/2, '#<_body'/3).
:- tool((>)/2, '#>_body'/3).
:- tool((>=)/2, '#>=_body'/3).
:- tool((=<)/2, '#=<_body'/3).

:- tool((=\=)/3, '#\\=_reif_body'/4).
:- tool((=:=)/3, '#=_reif_body'/4).
:- tool((<)/3, '#<_reif_body'/4).
:- tool((>)/3, '#>_reif_body'/4).
:- tool((>=)/3, '#>=_reif_body'/4).
:- tool((=<)/3, '#=<_reif_body'/4).
*/

:- tool('::'/2, '::_body'/3).
:- tool('#::'/2, '::_body'/3).

:- local reference(prob_handle,0).
:- local store(stats).

:- export gfd_var_print/2.
:- export gfd_copy_var/2.

load_gfd_solver(Arch) :- 
        get_flag(object_suffix, O),
        ( Arch = "x86_64_linux" ->
            concat_string(["gfd.", O], SolverObj),
            getcwd(Current),
            cd(Arch),
            block((load(SolverObj) -> cd(Current) ; cd(Current), fail), Tag,
                  (cd(Current), exit_block(Tag)))
        ;
            concat_string([Arch,/,"gfd.", O], SolverObj),
            load(SolverObj)
        ).

:- 
        get_flag(hostarch, Arch),
        load_gfd_solver(Arch),
        external(g_init/0, p_g_init),
        external(g_init_space_handle_c/1, p_g_init_space_handle_c),
        external(g_state_is_stable/1, p_g_state_is_stable),
        external(g_check_handle/3, p_g_check_handle),
        external(g_trail_undo_for_event/1, p_g_trail_undo_for_event),
        external(g_delete/1, p_g_delete),
        external(g_add_newbool/3, p_g_add_newbool),
        external(g_add_newvars_interval/4, p_g_add_newvars_interval),
        external(g_add_newvars_dom/3, p_g_add_newvars_dom),
        external(g_add_newvars_dom_handle/3, p_g_add_newvars_dom_handle),
        external(g_add_newvars_dom_union/4, p_g_add_newvars_dom_union),
        external(g_add_newvar_copy/3, p_g_add_newvar_copy),
        external(g_add_newvars_as_bool/3, p_g_add_newvars_as_bool),
        external(g_link_newbools/2, p_g_link_newbools),
        external(g_post_interval/4, p_g_post_interval),
        external(g_post_var_interval_reif/5, p_g_post_var_interval_reif),
        external(g_post_dom/3, p_g_post_dom),
        external(g_post_dom_handle/3, p_g_post_dom_handle),
        external(g_post_var_dom_reif/4, p_g_post_var_dom_reif),
        external(g_post_var_val_reif/4, p_g_post_var_val_reif),
        external(g_post_exclude_val/3, p_g_post_exclude_val),
        external(g_post_exclude_range/4, p_g_post_exclude_range),
        external(g_post_exclude_dom/3, p_g_post_exclude_dom),
        external(g_post_exclude_dom_handle/3, p_g_post_exclude_dom_handle),
        external(g_post_exclude_var_val/3, p_g_post_exclude_var_val),
        external(g_post_setvar/3, p_g_post_setvar),
        external(g_post_intrel_cstr/3, p_g_post_intrel_cstr),
        external(g_post_bool_connectives/3, p_g_post_bool_connectives),
        external(g_post_alldiff/3, p_g_post_alldiff),
        external(g_post_alldiff_offsets/4, p_g_post_alldiff_offsets),
        external(g_post_nvalues/4, p_g_post_nvalues),
        external(g_post_count/6, p_g_post_count),
        external(g_post_among/5, p_g_post_among),
        external(g_post_count_matches/5, p_g_post_count_matches),
        external(g_post_gcc/5, p_g_post_gcc),
        external(g_post_element/5, p_g_post_element),
        external(g_post_sorted2/4, p_g_post_sorted2),
        external(g_post_sorted/5, p_g_post_sorted),
        external(g_post_sequence/7, p_g_post_sequence),
        external(g_post_sequence_01/6, p_g_post_sequence_01),
        external(g_post_circuit/4, p_g_post_circuit),
        external(g_post_circuit_cost/7, p_g_post_circuit_cost),
        external(g_post_ham_path/6, p_g_post_ham_path),
        external(g_post_ham_path_cost/9, p_g_post_ham_path_cost),
        external(g_post_disj/4, p_g_post_disj),
        external(g_post_disjflex/5, p_g_post_disjflex),
        external(g_post_disjoint2/6, p_g_post_disjoint2),
        external(g_post_disjointflex2/8, p_g_post_disjointflex2),
        external(g_post_cumulative/6, p_g_post_cumulative),
        external(g_post_cumulativeflex/7, p_g_post_cumulativeflex),
        external(g_post_cumulatives/8, p_g_post_cumulatives),
        external(g_post_sum/5, p_g_post_sum),
        external(g_post_lin/6, p_g_post_lin),
        external(g_post_sum_reif/6, p_g_post_sum_reif),
        external(g_post_lin_reif/7, p_g_post_lin_reif),
        external(g_post_maxlist/4, p_g_post_maxlist),
        external(g_post_minlist/4, p_g_post_minlist),
        external(g_post_sqrt/4, p_g_post_sqrt),
        external(g_post_sq/4, p_g_post_sq),
        external(g_post_abs/4, p_g_post_abs),
        external(g_post_div/4, p_g_post_div),
        external(g_post_mult/5, p_g_post_mult),
        external(g_post_mod/4, p_g_post_mod),
        external(g_post_min2/5, p_g_post_min2),
        external(g_post_max2/5, p_g_post_max2),
        external(g_post_divmod/5, p_g_post_divmod),
        external(g_post_boolchannel/5, p_g_post_boolchannel),
        external(g_post_inverse/4, p_g_post_inverse),
        external(g_post_inverse_offset/6, p_g_post_inverse_offset),
        external(g_post_ordered/4, p_g_post_ordered),
        external(g_post_rel/5, p_g_post_rel),
        external(g_post_collection_rel/5, p_g_post_collection_rel),
        external(g_post_lex_order/5, p_g_post_lex_order),
        external(g_post_bin_packing/4, p_g_post_bin_packing),
        external(g_post_lwb/3, p_g_post_lwb),
        external(g_post_upb/3, p_g_post_upb),
        external(g_post_precede/4, p_g_post_precede),
        external(g_post_precede_chain/3, p_g_post_precede_chain),
        external(g_post_mem/3, p_g_post_mem),
        external(g_post_mem_reif/4, p_g_post_mem_reif),
        external(g_post_table/6, p_g_post_table),
        external(g_post_extensional/4, p_g_post_extensional),
        external(g_create_tupleset_handle/3, p_g_create_tupleset_handle),
        external(g_create_regdfa_handle/2, p_g_create_regdfa_handle),
        external(g_create_dfa_handle/4, p_g_create_dfa_handle),
        external(g_propagate/4, p_g_propagate),
        external(g_check_val_is_in_var_domain/3, p_g_check_val_is_in_var_domain),
        external(g_get_var_bounds/4, p_g_get_var_bounds),
        external(g_get_var_value/3, p_g_get_var_value),
        external(g_get_var_domain/3, p_g_get_var_domain),
        external(g_get_var_lwb/3, p_g_get_var_lwb),
        external(g_update_and_get_var_bound/5, p_g_update_and_get_var_bound),
        external(g_get_var_upb/3, p_g_get_var_upb),
        external(g_get_var_domain_size/3, p_g_get_var_domain_size),
        external(g_get_var_domain_width/3, p_g_get_var_domain_width),
        external(g_get_var_degree/3, p_g_get_var_degree),
        external(g_get_var_median/3, p_g_get_var_median),
        external(g_get_var_afc/3, p_g_get_var_afc),
        external(g_get_var_regret_lwb/3, p_g_get_var_regret_lwb),
        external(g_get_var_regret_upb/3, p_g_get_var_regret_upb),
        external(g_setup_search/10, p_g_setup_search),
        external(g_do_search/7, p_g_do_search),
        external(g_get_gfd_maxint/1, p_g_get_gfd_maxint),
        external(g_get_gfd_minint/1, p_g_get_gfd_minint),
        external(g_get_var_domain_handle/3, p_g_get_var_domain_handle),
        external(g_propagate_recompute/1, p_g_propagate_recompute),
        external(g_stop_caching/1, p_g_stop_caching),
        external(g_start_caching/1, p_g_start_caching),
        external(g_create_idxs_handle/2, p_g_create_idxs_handle),
        external(g_select/4, p_g_select),

	external(g_gecode_version/1, p_g_gecode_version),

        g_init,
	g_gecode_version(Version),
        printf(log_output, "Loaded Gecode solver %s%n", [Version]).

:- export struct(
        gfd_prob(
             cp_stamp,
             nvars,
             nevents,
             vars,
             prop,
             last_anc,
             space,
             events,
             events_tail
        )
   ).

:- export portray(gfd_prob/(property(arity) of gfd_prob), gfd_handle_tr_out/2, []).
:- export gfd_handle_tr_out/2.
gfd_handle_tr_out(gfd_prob{nvars:N},gfd_prob(nvars(N))).


/* gfd_space{} is the ECLiPSe level representation of a Gecode 
   computation state.  At the ECLiPSe level, a space is cloned 
   (producing a new gfd_space{}) when the cloning condition is 
   satisfied, and the new gfd_space{} is used for continuing
   execution, with the old gfd_space{} becoming an "ancestor" of
   the current one. The actual Gecode computation space pointed
   to by a gfd_space{} can change, e.g. when backtracking undo 
   changes to the computation state. In such cases, a new Gecode
   computation space is needed, and this is done by cloning the
   state from the nearest ancestor, and then recomputing the rest
*/
:- export struct( gfd_space(handle, stamp) ).

:- export struct(
        gfd(
            idx,
            bool,
            prob,
            any,
            set
        )
   ).

:- meta_attribute(gfd, [
        set_bounds:gfd_set_var_bounds/3,  
        get_bounds:gfd_get_var_bounds/3,  
        print:gfd_var_print/2,
        unify:gfd_unify/2,
        copy_term:gfd_copy_var/2,
	suspension_lists:[
	    any:(any of gfd),
	    min:(any of gfd),	% approximate alias for generic code
	    max:(any of gfd)	% approximate alias for generic code
	]
    ]).


:- local struct(options(interval_min,
                        interval_max,
                        array_size,
                        events_max,
                        cloning_distance)
               ).

valid_option_field(interval_min, interval_min of options).
valid_option_field(interval_max, interval_max of options).
valid_option_field(array_size,   array_size of options).
valid_option_field(cloning_distance, cloning_distance of options).
valid_option_field(events_max, events_max of options).

valid_option_value(interval_min, Min) :- integer(Min), Min >= gfd_minint.
valid_option_value(interval_max, Max) :- integer(Max), Max =< gfd_maxint.
valid_option_value(array_size, Size) :- integer(Size), Size > 0.
valid_option_value(cloning_distance, D) :- integer(D), D > 0.
valid_option_value(events_max, Max) :- integer(Max), Max > 0.

default_options(options{interval_min: -10000000,
                        interval_max:  10000000,
                        array_size: 100,
                        events_max: 2000,
                        cloning_distance: 2}
               ).


gfd_get_default(Option, Value) :-
        valid_option_field(Option, OptPos),
        get_options([], Defaults),
        arg(OptPos, Defaults, Value).

gfd_set_default(Option, Value) :-
        set_default_option(Option, Value),
        ( Option == cloning_distance ->
            % treat defaults that are acessed frequently specially
            % so that they can be accessed quickly -- accessing
            % a compiled term is faster than accessing a local variable
%            setval(cloning_distance, Value)
            compile_term(cloning_distance(Value))

        ; Option == interval_min ->
            gfd_get_default(interval_max, Max),
            compile_term(gfd_default_interval(Value,Max))

        ; Option == interval_max ->
            gfd_get_default(interval_min, Min),
            compile_term(gfd_default_interval(Min,Value))
        ; Option == events_max ->
            gfd_get_default(events_max, EMax),
            compile_term(events_max(EMax))
        ;                
            true
        ).


?- gfd_get_default(cloning_distance, Dist),
   compile_term(cloning_distance(Dist)),
   gfd_get_default(events_max, EMax),
   compile_term(events_max(EMax)),
   gfd_get_default(interval_min, Min),
   gfd_get_default(interval_max, Max),
   compile_term(gfd_default_interval(Min, Max)).

/*
gfd_default_interval(Min, Max) :-
        valid_option_field(interval_min, MinPos),
        valid_option_field(interval_max, MaxPos),
        get_options([], Defaults),
        arg(MinPos, Defaults, Min),
        arg(MaxPos, Defaults, Max).
*/

boolean_expr(E) :-
        reifiable_constraint(E), !.
boolean_expr(E) :-
        connective(E, _, _).


reifiable_constraint(sum(_Vs,_Rel,_S)) ?- !.
reifiable_constraint(scalar_product(_Cs,_Vs,_Rel,_P)) ?- !.
reifiable_constraint(mem(_Vs,_M)) ?- !.
/*reifiable_constraint(E) :-
        relation_constraint(E, _, _, _), !.
*/
reifiable_constraint((_ #= _)) ?- !.
reifiable_constraint((_ #\= _)) ?- !.
reifiable_constraint((_ #=< _)) ?- !.
reifiable_constraint((_ #< _)) ?- !.
reifiable_constraint((_ #>= _)) ?- !.
reifiable_constraint((_ #> _)) ?- !.
/*reifiable_constraint(E) :-
        domain_constraint(E,_, _), !.
*/
reifiable_constraint((_ :: _)) ?- !.
reifiable_constraint((_ #:: _)) ?- !.


domain_constraint((Vs0 :: Dom0), Vs, Dom) ?- !,
        Vs = Vs0, Dom = Dom0.
domain_constraint((Vs0 #:: Dom0), Vs, Dom) ?- !,
        Vs = Vs0, Dom = Dom0.

relation_constraint(X0 #= Y0, RelOp, X, Y) ?- !,
        RelOp = (#=), 
        X0 = X,
        Y0 = Y.
relation_constraint(X0 #\= Y0, RelOp, X, Y) ?- !,
        RelOp = (#\=), 
        X0 = X,
        Y0 = Y.
relation_constraint(X0 #< Y0, RelOp, X, Y) ?- !,
        RelOp = (#<), 
        X0 = X,
        Y0 = Y.
relation_constraint(X0 #=< Y0, RelOp, X, Y) ?- !,
        RelOp = (#=<), 
        X0 = X,
        Y0 = Y.
relation_constraint(X0 #> Y0, RelOp, X, Y) ?- !,
        RelOp = (#>), 
        X0 = X,
        Y0 = Y.
relation_constraint(X0 #>= Y0, RelOp, X, Y) ?- !,
        RelOp = (#>=), 
        X0 = X,
        Y0 = Y.


connective(X and Y, ConOp, Args) ?- !,
        ConOp = and,
        Args = [X,Y].
connective(X or Y, ConOp, Args) ?- !,
        ConOp = or,
        Args = [X,Y].
connective(X xor Y, ConOp, Args) ?- !,
        ConOp = xor,
        Args =[X,Y].
connective(X <=> Y, ConOp, Args) ?- !,
        ConOp = (<=>),
        Args = [X,Y].
connective(X => Y, ConOp, Args) ?- !,
        ConOp = (=>),
        Args = [X,Y].
connective(neg(X), ConOp, Args) ?- 
        ConOp = neg,
        Args = [X].

rel_op('#=').
rel_op('#\\=').
rel_op('#<').
rel_op('#=<').
rel_op('#>').
rel_op('#>=').

pl_rel_op('#=','=').
pl_rel_op('#\\=','\\=').
pl_rel_op('#<','<').
pl_rel_op('#=<','=<').
pl_rel_op('#>','>').
pl_rel_op('#>=','>=').

is_valid_rel_op(Rel0, Rel) :-
        atomic(Rel0), 
        ( rel_op(Rel0) -> Rel0 = Rel
        ; pl_rel_op(Rel, Rel0) -> true
        ; set_bip_error(6)
        ).
 
inline_op(X+Y, Out) ?- !, Out = X+Y.
inline_op(X-Y, Out) ?- !, Out = X-Y.
inline_op(X*Y, Out) ?- !, Out = X * Y.
inline_op(-X, Out) ?- !, Out = -X.
inline_op(isqrt(X), Out) ?- !, Out = isqrt(X).
inline_op(X^2, Out) ?- !, Out = sqr(X).
inline_op(sqr(X), Out) ?- !, Out = sqr(X).
inline_op(abs(X), Out) ?- !, Out = abs(X).
inline_op(X//Y, Out) ?- !, Out = X//Y.
inline_op(X rem Y, Out) ?- !, Out = X rem Y.
inline_op(max(X,Y), Out) ?- !, Out = max(X,Y).
inline_op(min(X,Y), Out) ?- !, Out = min(X,Y).
/* inline_op(element(Index, Collection, Out) -- not supported due to
   indexing differences
*/

/* aux_op are 'operators' that cannot be inlined in a Gecode expression.
   They are 'functional constraints' -- constraints written in a functional
   notation in an expression, without the last argument of the constraint,
   which is the value of the function. Only constraints that have a last
   argument as a domain variable can be functional constraint.
*/
aux_op(sum(_Vs,_Rel), Aux, S, _GS, Type, _ConLev) ?- !,
	Aux = linsum,
	Type = aux_cstr(S).
aux_op(scalar_product(_Cs,_Vs,_Rel), Aux, S, _GS, Type, _ConLev) ?- !,
	Aux = scalar_product,
	Type = aux_cstr(S).
aux_op(occurrences(_Value,_Vs), Aux, N, _GN, Type, _ConLev) ?- !,
	Aux = occurrences,
	Type = aux_cstr(N).
aux_op(count(_Value,_Vs,_Rel), Aux, N, _GN, Type, _ConLev) ?- !,
	Aux = count,
	Type = aux_cstr(N).
aux_op(among(_Values,_Vs,_Rel), Aux, N, _GN, Type, _ConLev) ?- !,
	Aux = among,
	Type = aux_cstr(N).
aux_op(count_matches(_Values,_Vs,_Rel), Aux, N, _GN, Type, _ConLev) ?- !,
	Aux = count_matches,
	Type = aux_cstr(N).
aux_op(element(_Idx,_Vs), Aux, Val, _GVal, Type, _ConLev) ?- !,
	Aux = element,
        Type = aux_cstr(Val).
aux_op(element_g(_Idx,_Vs), Aux, Val, _GVal, Type, _ConLev) ?- !,
	Aux = element_g,
        Type = aux_cstr(Val).
aux_op(mem(_Vs), Aux, Mem, _GVal, Type, _ConLev) ?- !,
	Aux = mem,
        Type = aux_cstr(Mem).
aux_op(circuit(_Succ, _CostMat), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = circuit,
        Type = aux_cstr(Cost).
aux_op(circuit(_Succ, _CostMat,_ArcCosts), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = circuit_arc,
        Type = aux_cstr(Cost).
aux_op(circuit_g(_Succ, _CostMat), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = circuit_g,
        Type = aux_cstr(Cost).
aux_op(circuit_g(_Succ, _CostMat,_ArcCosts), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = circuit_arc_g,
        Type = aux_cstr(Cost).
aux_op(circuit_offset(_Succ, _Off, _CostMat), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = circuitoff,
        Type = aux_cstr(Cost).
aux_op(circuit_offset(_Succ, _Off, _CostMat,_ArcCosts), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = circuitoff_arc,
        Type = aux_cstr(Cost).
aux_op(circuit_offset_g(_Succ, _Off, _CostMat), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = circuitoff_g,
        Type = aux_cstr(Cost).
aux_op(circuit_offset_g(_Succ, _Off, _CostMat,_ArcCosts), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = circuitoff_arc_g,
        Type = aux_cstr(Cost).
aux_op(ham_path(_S, _E, _Succ, _CostMat), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = ham_path,
        Type = aux_cstr(Cost).
aux_op(ham_path(_S, _E, _Succ, _CostMat,_ArcCosts), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = ham_path_arc,
        Type = aux_cstr(Cost).
aux_op(ham_path_g(_S, _E, _Succ, _CostMat), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = ham_path_g,
        Type = aux_cstr(Cost).
aux_op(ham_path_g(_S, _E, _Succ, _CostMat,_ArcCosts), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = ham_path_arc_g,
        Type = aux_cstr(Cost).
aux_op(ham_path_offset(_S, _E, _Succ, _Off, _CostMat), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = ham_pathoff,
        Type = aux_cstr(Cost).
aux_op(ham_path_offset(_S, _E, _Succ, _Off, _CostMat,_ArcCosts), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = ham_pathoff_arc,
        Type = aux_cstr(Cost).
aux_op(ham_path_offset_g(_S, _E, _Succ, _Off, _CostMat), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = ham_pathoff_g,
        Type = aux_cstr(Cost).
aux_op(ham_path_offset_g(_S, _E, _Succ, _Off, _CostMat,_ArcCosts), Aux, Cost, _GVal, Type, _ConLev) ?- !,
        Aux = ham_pathoff_arc_g,
        Type = aux_cstr(Cost).

%------------------------------------------------------------------------
% Support

% low level representation of variable
gfdvar(I,B,'_ivar'(I,B)).

% get the low-level representation of an existing gecide variable
get_gecode_var(_{gfd:Attr}, GV) ?-
        nonvar(Attr),
        Attr = gfd{idx:I,bool:B},
        gfdvar(I, B, GV).

get_gecode_var_or_int(_{gfd:Attr}, GV) ?-
        nonvar(Attr), !,
        Attr = gfd{idx:I,bool:B},
        gfdvar(I, B, GV).
get_gecode_var_or_int(I, GI) :-
        integer(I), !, 
        I = GI.
get_gecode_var_or_int(_, _) :-
        set_bip_error(6).

gfd_copy_var(_X{gfd:AttrX}, Copy) ?-
        ( var(AttrX) ->
            true
        ;
            AttrX = gfd{prob:H,idx:XIdx},
            restore_space_if_needed(H, _),
            H = gfd_prob{nvars:N0},
            new_gfdvar(Copy0, H, N0,N, _GCopy0),
            setarg(nvars of gfd_prob, H, N),
            post_new_event(copyvar(N,XIdx), H),
            Copy0 = Copy
        ).


gfd_unify(_Term, Attr) :-
        var(Attr).
gfd_unify(Term, Attr) :-
        compound(Attr),
        unify_term_gfd(Term, Attr).

:- mode unify_term_gfd(?, +).
unify_term_gfd(Y{gfd:AttrY}, AttrX) ?-
        unify_gfd_gfd(Y, AttrX, AttrY).
unify_term_gfd(X, Attr) :-
        integer(X),
        Attr = gfd{prob:H,set:S,idx:I},
        schedule_suspensions(any of gfd, Attr),
        ( S == [] -> 
            true 
        ; 
            post_new_event_no_wake(setvar(I, X), H)
        ).

unify_gfd_gfd(_Y, AttrX, AttrY) :-
        var(AttrY),
        AttrX = AttrY.
unify_gfd_gfd(_Y, AttrX, AttrY) :-
        nonvar(AttrY),
        AttrX = gfd{idx:IdxX,bool:BX,prob:H},
        AttrY = gfd{idx:IdxY,bool:BY,prob:H},
        % the variables must belong to the same problem, else fail
        ( IdxX == IdxY ->
            true   % same variable, do nothing
        ;
            % unless both variables have their own boolean already, make
            % sure they share the same boolean
            ( integer(BX), integer(BY) -> 
                true
            ;
                BX = BY
            ),
            merge_suspension_lists(any of gfd, AttrY, any of gfd, AttrX),
            % post an equality constraint for the two variables to gecode
            gfdvar(IdxX,BX,GX),
            gfdvar(IdxY,BY,GY),
            post_new_event_no_wake(post_rel(default,GX, (#=), GY), H)
        ).

gfd_set_var_bounds(_{gfd:Attr}, Lo0, Hi0) ?-
        nonvar(Attr),
        Lo is fix(Lo0),
        Hi is fix(Hi0),
        Attr = gfd{prob:H, idx:I, bool:BI},
        gfdvar(I,BI,GV),
        % note: changed to call wake as other bounds handler do so
        post_new_event(post_interval([](GV), Lo, Hi), H).

gfd_get_var_bounds(_{gfd:Attr}, Lo, Hi) ?-
        nonvar(Attr),
        Attr = gfd{prob:H, idx:I},
        % make sure there is a valid space
        restore_space_if_needed(H, SpH),
        g_get_var_bounds(SpH, I, Lo, Hi).

get_prob_handle(H) :-
        getval(prob_handle, H0),
        ( H0 = gfd_prob{} ->
            H0 = H
        ;
            new_prob_handle(H),
            setval(prob_handle, H)
        ).

/* returns prob_handle with current nvars, a common pattern, so optimise
   it a bit here
*/
get_prob_handle_nvars(H, NV) :-
        getval(prob_handle, H0),
        ( H0 = gfd_prob{nvars:NV} ->
            H = H0
        ;
            NV = 0,
            new_prob_handle(H),
            setval(prob_handle, H)
        ).

/* addto_varray(+ProbHandle, ++Idx, ?V)
     Add new problem variable V with index Idx to the variable array of 
     ProbHandle, expanding the array if required
*/
addto_varray(H, Idx, V) :-
  H = gfd_prob{vars:VArr},
  ( Idx > arity(VArr) ->
      expand_and_copy_array(VArr, NewVArr),
      setarg(vars of gfd_prob, H, NewVArr),
      arg(Idx, NewVArr, V)
  ;
      arg(Idx, VArr, V)
  ).

expand_and_copy_array(Old, New) :-
        arity(Old, OldSize),
        NewSize is OldSize*2,
        dim(New, [NewSize]),
        ( foreacharg(A, Old, Idx), param(New) do
            arg(Idx, New, A)
        ).

create_and_add_default_gfdvar(V, H) :-
        H = gfd_prob{nvars:N0},
        gfd_default_interval(Min, Max),
        new_gfdvar(V, H, N0,N1, _),
        do_update_newvars_with_domain_interval(H, N1, Min, Max).

is_not_boolvar(_{gfd:Attr}) ?-
        nonvar(Attr), !,
        Attr = gfd{bool:Link},
        var(Link).
is_not_boolvar(_).

link_var_to_boolvar(V, H) :-
        ( (var(V), is_not_boolvar(V)) ->
	    get_gecode_attr(V, Attr),
            Attr = gfd{idx:Idx,bool:BIdx},
            % never done on its own, so no need to wake
            post_new_event_no_wake(newbool(Idx,BIdx), H)
        ;
            true
        ).
        
update_vars_for_gecode(N0, N, Bs, H, Min, Max) :-
        ( N > N0 ->
            % have new variables, add them (and link any booleans to their var)
            do_update_newvars_with_domain_interval(H, N, Min, Max)
        ;
            true
        ),
        ( foreach(B, Bs), param(H) do % link the boolean vars
            link_var_to_boolvar(B, H) 
        ).


get_gecode_attr(_{gfd:Attr0}, Attr) ?-
        nonvar(Attr0),
        Attr = Attr0.

add_gecode_attr(X{gfd:Attr}, H, Idx, BI) ?-
        var(Attr),
        new_gecode_attr(X, H, Idx, BI, Attr).
add_gecode_attr(X, H, Idx, BI) :-
        free(X),
        new_gecode_attr(X, H, Idx, BI, _Attr).

get_gecode_domain(X, Domain) :-
        get_gecode_attr(X, Attr),
        Attr = gfd{idx:Idx,prob:gfd_prob{space:gfd_space{handle:SpH}}},
        g_get_var_domain(SpH, Idx, Domain).

gfd_var_print(X, Domain) :-
        get_gecode_domain(X, Domain).

:- mode new_gecode_attr(?,+,+,?,-).
new_gecode_attr(X, H, N, BN, Attr) :-
        Attr = gfd{prob:H,idx:N, bool:BN},
        init_suspension_list(any of gfd, Attr),
        add_attribute(X, Attr, gfd).


new_prob_handle(H) :-
        gfd_get_default(array_size, VSz),
        dim(VArr, [VSz]),
        H = gfd_prob{nvars:0,last_anc:[],space:Sp,vars:VArr,
                 nevents:0, prop:Susp},
        % setarg/3 instead of arg/3 is used to initialise the events list
        % to make sure the variable tail is not physically allocated
        % inside the gfd_prob structure (idea taken from notify_port.ecl)
        setarg(events_tail of gfd_prob, H, Tail),
        setarg(events of gfd_prob, H, Tail),
        timestamp_update(H, cp_stamp of gfd_prob),
        make_suspension(gfd_do_propagate(H), 9, Susp),
        new_space_handle(Sp).

new_space_handle(Sp) :-
        Sp = gfd_space{handle:SH},
        timestamp_init(Sp, stamp of gfd_space),
        g_init_space_handle_c(SH).


% create a new gfd variable at the ECLiPSe level. Note: variable need to
% be added to Gecode
new_gfdvar(V, H, N0, N, GV) :-
        N is N0 + 1,
        gfdvar(N, BN, GV),
        add_gecode_attr(V, H, N, BN),  % may fail!
        addto_varray(H, N, V).


% converts a compound term to a list: a list is left as a list
% (without flattening), other compound terms are converted into a list
% with its arguments as the elements of the list, i.e. the term is not
% flattened like in collection_to_list
check_compound_to_list(Cs, L) :-
        check_nonvar(Cs), 
        ( Cs = [_|_] ->
            Cs = L
        ; compound(Cs) ->
            Cs =.. [_|L]
        ;
            set_bip_error(5)  % type error
        ).

        
% type/range checking used in constraints

check_integer(I) :-
        (integer(I) -> true ; set_bip_error(5)).

check_nonnegative(I) :-   % assumes I is a number
        (I >= 0 -> true ; set_bip_error(6)).


check_atom(A) :-
        (atom(A) -> true ; set_bip_error(5)).

check_nonvar(A) :-
        (nonvar(A) -> true ; set_bip_error(4)). % instantiation fault

check_collection_to_list(C, L) :-
        (collection_to_list(C, L) -> true ; set_bip_error(5)).


%------------------------------------------------------------------------
% Expression support

post_connectives(Conn, ConLev, Module) :-
        get_prob_handle_nvars(H, N0),
        ec_to_gecode_bool_expr1(Conn, H, N0,N, [],Bs, Auxs0,AuxsT, GConn, ConLev, Module),
%        ec_to_gecode_connectives1(Conn, H, N0,N, [],Bs, Auxs0,AuxsT, GConn, Module),
        !,
        gfd_default_interval(Min, Max),
        update_vars_for_gecode(N0, N, Bs, H, Min, Max),
        post_new_event_with_aux([post_bool_connectives(ConLev,GConn)|Auxs0],AuxsT, H).
post_connectives(Conn, _ConLev, _Module) :-
        get_bip_error(E),
        error(E, Conn).


:- tool((and)/2, and_body/3).
:- tool((or)/2, or_body/3).
:- tool((xor)/2, xor_body/3).
:- tool(neg/1, neg_body/2).
:- tool('<=>'/2, '<=>_body'/3).
:- tool('=>'/2, '=>_body'/3).

:- tool((and)/3, and_reif_body/4).
:- tool((or)/3, or_reif_body/4).
:- tool((xor)/3, xor_reif_body/4).
:- tool(neg/2, neg_reif_body/3).
:- tool('<=>'/3, '<=>_reif_body'/4).
:- tool('=>'/3, '=>_reif_body'/4).

:- tool(and_c/3, and_c/4).
:- tool(or_c/3, or_c/4).
:- tool(xor_c/3, xor_c/4).
:- tool(neg_c/2, neg_c/3).
:- tool('<=>_c'/3, '<=>_c'/4).
:- tool('=>_c'/3, '=>_c'/4).

:- tool(and_reif_c/4, and_reif_c/5).
:- tool(or_reif_c/4, or_reif_c/5).
:- tool(xor_reif_c/4, xor_reif_c/5).
:- tool(neg_reif_c/3, neg_reif_c/4).
:- tool('<=>_reif_c'/4, '<=>_reif_c'/5).
:- tool('=>_reif_c'/4, '=>_reif_c'/5).

:- tool(among/4, among_body/5).
:- tool(among_c/5, among_c/6).


and_body(EX, EY, Module) :-
        and_c(EX, EY, default, Module).

and_c(EX, EY, ConLev, Module) :-
        post_connectives((EX and EY), ConLev, Module).

or_body(EX, EY, Module) :-
        or_c(EX, EY, default, Module).

or_c(EX, EY, ConLev, Module) :-
        post_connectives((EX or EY), ConLev, Module).

xor_body(EX, EY, Module) :-
        xor_c(EX, EY, default, Module).

xor_c(EX, EY, ConLev, Module) :-
        post_connectives((EX xor EY), ConLev, Module).

'<=>_body'(EX, EY, Module) :-
        '<=>_c'(EX, EY, default, Module).

'<=>_c'(EX, EY, ConLev, Module) :-
        post_connectives((EX <=> EY), ConLev, Module).

'=>_body'(EX, EY, Module) :-
        '=>_c'(EX, EY, default, Module).

'=>_c'(EX, EY, ConLev, Module) :-
        post_connectives((EX => EY), ConLev, Module).


neg_body(EX, Module) :-
        neg_c(EX, default, Module).

neg_c(EX, ConLev, Module) :-
        post_connectives(neg(EX), ConLev, Module).


and_reif_body(EX, EY, Bool, Module) :-
        and_reif_c(EX, EY, Bool, default, Module).

and_reif_c(EX, EY, Bool, ConLev, Module) :-
        post_connectives((Bool <=> (EX and EY)), ConLev, Module).

or_reif_body(EX, EY, Bool, Module) :-
        or_reif_c(EX, EY, Bool, default, Module).

or_reif_c(EX, EY, Bool, ConLev, Module) :-
        post_connectives((Bool <=> (EX or EY)), ConLev, Module).

xor_reif_body(EX, EY, Bool, Module) :-
        xor_reif_c(EX, EY, Bool, default, Module).

xor_reif_c(EX, EY, Bool, ConLev, Module) :-
        post_connectives((Bool <=> (EX xor EY)), ConLev, Module).

'<=>_reif_body'(EX, EY, Bool, Module) :-
        '<=>_reif_c'(EX, EY, Bool, default, Module).

'<=>_reif_c'(EX, EY, Bool, ConLev, Module) :-
        post_connectives((Bool <=> (EX <=> EY)), ConLev, Module).

'=>_reif_body'(EX, EY, Bool, Module) :-
        '=>_reif_c'(EX, EY, Bool, default, Module).

'=>_reif_c'(EX, EY, Bool, ConLev, Module) :-
        post_connectives((Bool <=> (EX => EY)), ConLev, Module).

neg_reif_body(EX, Bool, Module) :-
        neg_reif_c(EX, Bool, default, Module).

neg_reif_c(EX, Bool, ConLev, Module) :-
        post_connectives((Bool <=> neg(EX)), ConLev, Module).


'#\\=_body'(EX, EY, Module) :-
        '#\\=_c'(EX, EY, default, Module).

'#\\=_c'(EX, EY, ConLev, Module) :-
        post_rel_cstr((#\=), EX, EY, ConLev, Module).

'#=_body'(EX, EY, Module) :-
        '#=_c'(EX, EY, default, Module).

'#=_c'(EX, EY, ConLev, Module) :-
        % optimisation for top-level aux. expressions
        ( (aux_op(EX, EXTemp, Res, GRes, EXType, ConLev),
	   aux_op(EY, EYTemp, Res, GRes, EYType,ConLev)) ->
            % EX and EY both aux 
	    get_prob_handle_nvars(H, N0),
	    new_gfdvar(Res, H, N0,N1, GRes),
	    ec_to_gecode_aux_op1(EXType, EX, H, N1,N2, [],Bs1, Auxs0,Auxs1, EXTemp, ConLev, Module),
	    ec_to_gecode_aux_op1(EYType, EY, H, N2,N3, Bs1,Bs, Auxs1,AuxsT, EYTemp, ConLev, Module),
            !,
            gfd_default_interval(Min, Max),
	    update_vars_for_gecode(N0, N3, Bs, H, Min, Max),
            post_new_event_with_aux(Auxs0,AuxsT, H)

        ;
            !, post_rel_cstr((#=), EX, EY, ConLev, Module)
        ).
'#=_c'(EX, EY, _ConLev, _Module) :-
        get_bip_error(E),
        error(E, (EX #= EY)).


'#<_body'(EX, EY, Module) :-
        '#<_c'(EX, EY, default, Module).

'#<_c'(EX, EY, ConLev, Module) :-
        post_rel_cstr((#<), EX, EY, ConLev, Module).

'#>_body'(EX, EY, Module) :-
        '#>_c'(EX, EY, default, Module).

'#>_c'(EX, EY, ConLev, Module) :-
        post_rel_cstr((#>), EX, EY, ConLev, Module).

'#>=_body'(EX, EY, Module) :-
        '#>=_c'(EX, EY, default, Module).

'#>=_c'(EX, EY, ConLev, Module) :-
        post_rel_cstr((#>=), EX, EY, ConLev, Module).

'#=<_body'(EX, EY, Module) :-
        '#=<_c'(EX, EY, default, Module).

'#=<_c'(EX, EY, ConLev, Module) :-
        post_rel_cstr((#=<), EX, EY, ConLev, Module).

'#\\=_reif_body'(EX, EY, Bool, Module) :-
        '<=>_c'((EX #\= EY), Bool, default, Module).

'#=_reif_body'(EX, EY, Bool, Module) :-
        '<=>_c'((EX #= EY), Bool, default, Module).

'#<_reif_body'(EX, EY, Bool, Module) :-
        '<=>_c'((EX #< EY), Bool, default, Module).

'#>_reif_body'(EX, EY, Bool, Module) :-
        '<=>_c'((EX #> EY), Bool, default, Module).

'#>=_reif_body'(EX, EY, Bool, Module) :-
        '<=>_c'((EX #>= EY), Bool, default, Module).

'#=<_reif_body'(EX, EY, Bool, Module) :-
        '<=>_c'((EX #=< EY), Bool, default, Module).

'#\\=_reif_c'(EX, EY, Bool, ConLev, Module) :-
        '<=>_c'((EX #\= EY), Bool, ConLev, Module).

'#=_reif_c'(EX, EY, Bool, ConLev, Module) :-
        '<=>_c'((EX #= EY), Bool, ConLev, Module).

'#<_reif_c'(EX, EY, Bool, ConLev, Module) :-
        '<=>_c'((EX #< EY), Bool, ConLev, Module).

'#>_reif_c'(EX, EY, Bool, ConLev, Module) :-
        '<=>_c'((EX #> EY), Bool, ConLev, Module).

'#>=_reif_c'(EX, EY, Bool, ConLev, Module) :-
        '<=>_c'((EX #>= EY), Bool, ConLev, Module).

'#=<_reif_c'(EX, EY, Bool, ConLev, Module) :-
        '<=>_c'((EX #=< EY), Bool, ConLev, Module).



post_rel_cstr(RelOp, EX, EY, ConLev, Module) :-
        get_prob_handle_nvars(H, N0),
        ec_to_gecode_expr1(EX, H, N0,N1, [],Bs1, Auxs0,Auxs1, GEX, ConLev, Module),
        ec_to_gecode_expr1(EY, H, N1,N2, Bs1,Bs2, Auxs1,Auxs2, GEY, ConLev, Module),
        construct_relcstr_event1(RelOp, EX, EY, GEX, GEY, H, N2,N, Bs2,Bs,
                                 Auxs2,AuxsT, Event, ConLev), 
        !,
        gfd_default_interval(Min, Max),
        update_vars_for_gecode(N0, N, Bs, H, Min, Max),
        post_new_event_with_aux([Event|Auxs0],AuxsT, H).
post_rel_cstr(RelOp, EX, EY, _ConLev, _Module) :-
        get_bip_error(E),
        Goal =.. [RelOp, EX, EY],
        error(E, Goal).


construct_relcstr_event1((#=), EX, EY, GEX, GEY, H, N0,N, Bs0,Bs, Auxs0,AuxsT, Event, ConLev) ?-
        !,
        ( boolean_expr(EX) ->
            ( boolean_expr(EY) ->
                Bs = Bs0,
                N = N0,
                Event = post_bool_connectives(ConLev,GEX<=>GEY),
                Auxs0 = AuxsT
            ; var(EY) ->
                Bs = [EY|Bs0],
                N = N0,
                Auxs0 = AuxsT,
                Event = post_bool_connectives(ConLev,GEY<=>GEX)
            ; integer(EY) ->
                Bs = Bs0,
                N = N0,
                Auxs0 = AuxsT,
                ( EY == 1 ->
                    Event = post_bool_connectives(ConLev,GEX)
                ; EY == 0 ->
                    % let Gecode do the optimisation if EY == 0
                    Event = post_bool_connectives(ConLev,GEY<=>GEX)
                ;
                    fail   % not boolean
                )
            ; 
                new_gfdvar(BVar, H, N0,N, GBVar),
                Bs = [BVar|Bs0],
                Event = post_bool_connectives(ConLev,GBVar<=>GEX),
                Auxs0 = [post_rc(ConLev,GBVar #= GEY)|AuxsT]
            )
        ; boolean_expr(EY) ->
            ( var(EX) ->
                Bs = [EX|Bs0],
                N = N0,
                Auxs0 = AuxsT,
                Event = post_bool_connectives(ConLev,GEX<=>GEY)
            ; integer(EX) ->
                Bs = Bs0,
                N = N0,
                Auxs0 = AuxsT,
                ( EX == 1 ->
                    Event = post_bool_connectives(ConLev,GEY)
                ; EX == 0 ->
                    Event = post_bool_connectives(ConLev,GEX<=>GEY)
                ;
                    fail  % not boolean
                )
            ;
                new_gfdvar(BVar, H, N0,N, GBVar),
                Bs = [BVar|Bs0],
                Event = post_bool_connectives(ConLev,GBVar<=>GEY),
                Auxs0 = [post_rc(ConLev,GBVar #= GEX)|AuxsT]
            )
        ;
            Bs = Bs0,
            N = N0,
            Auxs0 = AuxsT,
            Event = post_rc(ConLev,GEX #= GEY)
        ).
construct_relcstr_event1(RelOp, EX, EY, GEX, GEY, H, N0,N, Bs0,Bs, Auxs0,AuxsT, Event, ConLev) ?-
        ( boolean_expr(EX) ->
            new_gfdvar(BXV, H, N0,N1, GBXV),
            Bs1 = [BXV|Bs0],
            Auxs0 = [post_bool_connectives(ConLev,GBXV<=>GEX)|Auxs1],
            GNewX = GBXV
        ;
            N1 = N0,
            Bs1 = Bs0,
            Auxs0 = Auxs1,
            GNewX = GEX
        ),
        ( boolean_expr(EY) ->
            new_gfdvar(BYV, H, N1,N, GBYV),
            Bs = [BYV|Bs1],
            Auxs1 = [post_bool_connectives(ConLev,GBYV<=>GEY)|AuxsT],
            GNewY = GBYV
        ;
            N = N1,
            Bs = Bs1,
            Auxs1 = AuxsT,
            GNewY = GEY
        ),
        RC =.. [RelOp, GNewX, GNewY],
        Event = post_rc(ConLev,RC).


ec_to_gecode_var(V, H, GV) :-
        H = gfd_prob{nvars:N0},
        ec_to_gecode_var1(V, H, N0,N, GV),
        update_gecode_with_default_newvars(H, N0, N).

update_gecode_with_default_newvars(H, N0, N) :-
        ( N > N0 ->
            % V is a new gecode var, add it to gecode
            gfd_default_interval(Min, Max),
            do_update_newvars_with_domain_interval(H, N, Min, Max)
        ;
            true
        ).

ec_to_gecode_oldvarlist_bounds1(Vs, H, Min, Max, N0,N, OldGVs) :-
        ( foreach(V, Vs),
          fromto(N0, N1,N2, N),
          param(H, Min, Max),
          fromto(OldGVs, GVs1,GVs2, [])
        do
            ( integer(V) ->
                V >= Min,
                V =< Max,
                N1 = N2,
                GVs1 = GVs2
            ; var(V) ->
                ( get_gecode_attr(V, gfd{idx:I,bool:B}) ->
                    gfdvar(I, B, GV),
                    GVs1 = [GV|GVs2],
                    N1 = N2
                ; /* new gfd var */
                    new_gfdvar(V, H, N1,N2, _GV),
                    GVs1 = GVs2
                )
            ;
                fail
            )
        ).

ec_to_gecode_varlist(L, H, GL, HasVar) :-
        H = gfd_prob{nvars:N0},
        ec_to_gecode_varlist1(L, H, N0,N, GL, HasVar),
        update_gecode_with_default_newvars(H, N0, N).

ec_to_gecode_varlist1(L, H, N0,N, GL, HasVar) :-
        ( foreach(V, L), 
          fromto(N0, N1,N2, N),
          param(H, HasVar),
          foreach(GV, GL)
        do
            ( var(V) ->
                HasVar = 1, % HasVar = 1 if there is at least 1 var
                ec_to_gecode_var1(V, H, N1,N2, GV)
            ; integer(V) ->
                GV = V,
                N1 = N2
            ;
                printf(error, "Integer or gfd variable expected for: %w%n", [V]),
                set_bip_error(5)
            )
        ).


ec_to_gecode_multivarlists(Ls, H, GLs) :-
        H = gfd_prob{nvars:NV0},
        ec_to_gecode_multivarlists1(Ls, H, NV0, NV, GLs),
        update_gecode_with_default_newvars(H, NV0, NV).

ec_to_gecode_multivarlists1(Ls, H, NV0,NV, GLs) :-
        ( foreach(L, Ls),
          param(H),
          fromto(NV0, NV1,NV2, NV),
          foreach(GL, GLs)
        do
            ec_to_gecode_varlist1(L, H, NV1,NV2, GL, _)
        ).


update_gecode_with_boolvars(NewVs, NBVs, N0, N1, HasVar, H) :-
        ( nonvar(HasVar) ->
            ( NBVs \== [] ->
                % new bool vars for existing int vars
                NBVArr =.. [[]|NBVs],
                Es = [connectnewbools(NBVArr)|Es1]
            ;
                Es = Es1
            ),
            ( N1 > N0 ->
                % new vars
                setarg(nvars of gfd_prob, H, N1),
                NewVArr =.. [[]|NewVs],
                Es1 = [newboolvars(N1, NewVArr)|EsTail]
            ;
                Es1 = EsTail
            ),
            ( nonvar(Es) ->
                post_new_event_with_aux(Es, EsTail, H)
            ;
                true  % no new bool vars
            )
        ;
            true  % no vars
        ).

ec_to_gecode_boolvarlist1(L, H, NewVs, NBVs, N0,N, GL, HasVar) :-
        ( foreach(V, L),
          fromto(NewVs, NewVs0,NewVs1, []),
          fromto(NBVs, NBVs0,NBVs1, []),
          fromto(N0, N1,N2, N),
          param(H, HasVar, N0),
          foreach(GV, GL)
        do
            ( var(V) ->
                HasVar = 1,
                ec_to_gecode_var1(V, H, N1,N2, GV),
                ( N1 == N2 -> 
                    % not a new domain variable
                    % need to check if already linked to a bool var
                    NewVs0 = NewVs1,
                    gfdvar(VI,BI, GV),
                    (integer(BI) ->
                        % already linked to a bool var
                        NBVs0 = NBVs1
                    ; VI > N0 -> 
                        % is a just created domain variable 
                        NBVs0 = NBVs1
                    ;
                        NBVs0 = [GV|NBVs1]
                    )
                ;
                    NewVs0 = [GV|NewVs1],
                    NBVs0 = NBVs1
                )
            ; integer(V) ->
                V >= 0,
                V =< 1,  % fails if V (a boolean) is out of range
                GV = V,
                N1 = N2,
                NewVs0 = NewVs1,
                NBVs0 = NBVs1
            ;
                printf(error, "0/1 integer or gfd variable expected for: %w%n", [V]),
                set_bip_error(5)
            )
        ).


ec_to_gecode_arith_exprlist1(List, H, Inline, N0,N, Bs0,Bs, Auxs0,Auxs, GList, ConLev, Module) :-
        ( foreach(E,List), fromto(N0, N1,N2, N), 
          fromto(Auxs0, As1,As2, Auxs), fromto(Bs0, Bs1,Bs2, Bs),
          param(H,Module,Inline,ConLev),
          foreach(GE,GList)
        do
            ( var(E) ->
                As2 = As1,
                Bs2 = Bs1,
                ec_to_gecode_var1(E, H, N1,N2, GE)
            ; integer(E) ->
                N2 = N1,
                GE = E,
                Bs2 = Bs1,
                As2 = As1
	    ; E = subscript(T,S) ->
	        subscript(T, S, V),
		As2 = As1,
		Bs2 = Bs1,
		ec_to_gecode_var1(V, H, N1,N2, GE)
           
            ; Inline == 0 -> 
                % only parse for subexpressions if non-inline allowed
                new_gfdvar(_NewV, H, N1,N3, GE),
                ec_to_gecode_arith_expr1(E, H, Inline, N3,N2, Bs1,Bs2, As1,As3, GExpr, ConLev, Module),
                As3 = [post_rc(ConLev,GE #= GExpr)|As2]
            ;
                printf(error, "Non-inlined expressions not allowed: %w%n", [E]),
                set_bip_error(21)
            )
        ).



ec_to_gecode_aux_op1(args(Specs), AuxOp, H, N0,N, Bs0,Bs, Auxs0,AuxsT, EventTemp, ConLev, Module) ?- !,
        % operator has multiple arguments of different types, with the ones
        % that can be expressions specified in the list Specs 
        Auxs0 = [EventTemp|Auxs1],
        ( foreach(Spec, Specs),
	  fromto(N0, N1,N2, N), 
	  fromto(Bs0, Bs1,Bs2, Bs), 
	  fromto(Auxs1, Auxs2,Auxs3, AuxsT),
	  param([EventTemp,AuxOp,H,ConLev,Module])
	do 
            ( Spec = list(AuxArg,EvArg) ->
                arg(AuxArg, AuxOp, Collect),
                check_collection_to_list(Collect, List),
                ec_to_gecode_arith_exprlist1(List, H, 0, N1,N2, Bs1,Bs2, Auxs2,Auxs3, GList, ConLev, Module),
                GArray =.. [[]|GList],
                arg(EvArg, EventTemp, GArray)
            ; Spec = var(AuxArg,EvArg) ->
                arg(AuxArg, AuxOp, Arg),
                arg(EvArg, EventTemp, GArg),
                ( var(Arg) ->
                    Auxs3 = Auxs2,
                    Bs2 = Bs1,
                    ec_to_gecode_var1(Arg, H, N1,N2, GArg)
                ; integer(Arg) ->
                    Auxs3 = Auxs2,
                    Bs1 = Bs2,
                    N1 = N2,
                    GArg = Arg
                ; % Arg is an expression, process it and add an auxillary
                  % constraint linking it to the argument variable in EventTemp
                    new_gfdvar(_NewV, H, N1,N3, GArg),
                    Auxs2 = [post_rc(ConLev, GArg #= GExpr)|Auxs4],
                    ec_to_gecode_arith_expr1(Arg, H, 0, N3,N2, Bs1,Bs2,
                                             Auxs4,Auxs3, GExpr, ConLev, Module)
                )
            ;
                set_bip_error(21)
            )
	     
        ).
ec_to_gecode_aux_op1(aux_cstr(Val), Cstr, H, N0,N, Bs0,Bs, Auxs0,AuxsT, EventTemp, ConLev, Module) ?- !,
        cstr_aux_events1(EventTemp, Cstr, Val, H, N0,N, Bs0,Bs, Auxs0,
                         AuxsT, ConLev, Module).


cstr_aux_events1(linsum, sum(Vs,Rel), Sum, H, N0,N, Bs0,Bs, Auxs0, AuxsT, ConLev, _Module) :-
        Bs0 = Bs,
        linsum_event1(Vs, Rel, Sum, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(scalar_product, scalar_product(Cs,Vs,Rel), Prod, H, N0,N, Bs0,Bs, Auxs0, AuxsT, ConLev, _Module) :-
        Bs0 = Bs,
        scalar_product_event1(Cs, Vs, Rel, Prod, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(element, element(Idx,Vs), Val, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        element_events1(Idx, Vs, Val, ecl, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(element_g, element_g(Idx,Vs), Val, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        element_events1(Idx, Vs, Val, gc, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(occurrences, occurrences(Value,Vs), Count, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        count_events1(Value, Vs, '#=', Count, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(count, count(Value,Vs,Rel), Val, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        count_events1(Value, Vs, Rel, Val, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(among, among(Values,Vs,Rel), Val, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, Module) :-
	Bs0 = Bs,
        among_events1(Values, Vs, Rel, Val, ConLev, H, N0,N, Auxs0, AuxsT, Module).
cstr_aux_events1(count_matches, count_matches(Values,Vs,Rel), Val, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        count_matches_events1(Values, Vs, Rel, Val, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(mem, mem(Vs), Mem, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        mem_events1(Vs, Mem, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(circuit, circuit(Succ, CostMat), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        circuit_aux_events1(Succ, 1, CostMat, [], Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(circuit_g, circuit_g(Succ, CostMat), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        circuit_aux_events1(Succ, 0, CostMat, [], Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(circuit_arc, circuit(Succ, CostMat, ArcCs), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        circuit_aux_events1(Succ, 1, CostMat, ArcCs, Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(circuit_arc_g, circuit_g(Succ, CostMat, ArcCs), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        circuit_aux_events1(Succ, 0, CostMat, ArcCs, Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(circuitoff, circuit_offset(Succ, Off0, CostMat), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        Off is Off0 + 1,
        circuit_aux_events1(Succ, Off, CostMat, [], Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(circuitoff_g, circuit_offset_g(Succ, Off, CostMat), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        circuit_aux_events1(Succ, Off, CostMat, [], Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(circuitoff_arc, circuit(Succ, Off0, CostMat, ArcCs), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        Off is Off0 + 1,
        circuit_aux_events1(Succ, Off, CostMat, ArcCs, Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(circuitoff_arc_g, circuit_offset_g(Succ, Off, CostMat, ArcCs), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        circuit_aux_events1(Succ, Off, CostMat, ArcCs, Cost, ConLev, H, N0,N, Auxs0,AuxsT).

cstr_aux_events1(ham_path, ham_path(Start, End, Succ, CostMat), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        ham_path_aux_events1(Start, End, Succ, 1, CostMat, [], Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(ham_path_g, ham_path_g(Start, End, Succ, CostMat), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        ham_path_aux_events1(Start, End, Succ, 0, CostMat, [], Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(ham_path_arc, ham_path(Start, End, Succ, CostMat, ArcCs), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        ham_path_aux_events1(Start, End, Succ, 1, CostMat, ArcCs, Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(ham_path_arc_g, ham_path_g(Start, End, Succ, CostMat, ArcCs), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        ham_path_aux_events1(Start, End, Succ, 0, CostMat, ArcCs, Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(ham_pathoff, ham_path_offset(Start, End, Succ, Off0, CostMat), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        Off is Off0 + 1,
        ham_path_aux_events1(Start, End, Succ, Off, CostMat, [], Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(ham_pathoff_g, ham_path_offset_g(Start, End, Succ, Off, CostMat), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        ham_path_aux_events1(Start, End, Succ, Off, CostMat, [], Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(ham_pathoff_arc, ham_path(Start, End, Succ, Off0, CostMat, ArcCs), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        Off is Off0 + 1,
        ham_path_aux_events1(Start, End, Succ, Off, CostMat, ArcCs, Cost, ConLev, H, N0,N, Auxs0,AuxsT).
cstr_aux_events1(ham_pathoff_arc_g, ham_path_offset_g(Start, End, Succ, Off, CostMat, ArcCs), Cost, H, N0,N, Bs0,Bs, Auxs0,AuxsT, ConLev, _Module) :-
	Bs0 = Bs,
        ham_path_aux_events1(Start, End, Succ, Off, CostMat, ArcCs, Cost, ConLev, H, N0,N, Auxs0,AuxsT).



ec_to_gecode_reified1(scalar_product(Cs,Vs,Rel,S), H, N0,N, Bs0,Bs, Auxs0,AuxsT, Event, GBool, ConLev, _Module) :-
        !,
        Auxs0 = AuxsT,
        scalar_product_reif_event1(Cs,Vs, Rel, S, _Bool, ConLev, H, N0,N, Bs0,Bs, Event, GBool).
ec_to_gecode_reified1(sum(Vs,Rel,S), H, N0,N, Bs0,Bs, Auxs0,AuxsT, Event, GBool, ConLev, _Module) :-
        !,
        Auxs0 = AuxsT,
        linsum_reif_event1(Vs, Rel, S, _Bool, ConLev, H, N0,N, Bs0,Bs, Event, GBool).
ec_to_gecode_reified1(mem(Vs,M), H, N0,N, Bs0,Bs, Auxs0,AuxsT, Event, GBool, ConLev, _Module) :-
        !,
        Auxs0 = AuxsT,
        mem_reif_event1(Vs, M, _Bool, ConLev, H, N0,N, Bs0,Bs, Event, GBool).
ec_to_gecode_reified1(C, H, N0,N, Bs0,Bs, Auxs0,AuxsT, Event, GBool, _ConLev, Module) :-
        domain_constraint(C, Vs, Dom), !,
        Auxs0 = AuxsT,
        ec_to_gecode_domain_reified1(Vs, Dom, _Bool, H, N0,N, Bs0,Bs, Event, GBool, Module).
ec_to_gecode_reified1(C, H, N0,N, Bs0,Bs, Auxs0,AuxsT, Event, GBool, ConLev, Module) :-
        relation_constraint(C, RelOp, E1, E2), !,
        new_gfdvar(Bool, H, N0, N1, GBool),
        Bs1 = [Bool|Bs0],
        Event = post_bool_connectives(ConLev,GBool<=>GRC),
        ec_to_gecode_reifiedrc1(RelOp, E1, E2, H, N1,N, Bs1,Bs, Auxs0,AuxsT, GRC, ConLev, Module).


ec_to_gecode_reifiedrc1(RelOp, E1, E2, H, N0,N, Bs0,Bs, Auxs0,AuxsT, GRC, ConLev, Module) :-
        % the subexpressions must be inline here
        ec_to_gecode_arith_expr1(E1, H, 1, N0,N1, Bs0,Bs1, Auxs0,Auxs1, GE1, ConLev, Module),
        ec_to_gecode_arith_expr1(E2, H, 1, N1,N, Bs1,Bs, Auxs1,AuxsT, GE2, ConLev, Module),
        GRC =.. [RelOp,GE1,GE2].


ec_to_gecode_arith_expr(E, H, Auxs0,AuxsT, GE, ConLev, Module) :-
        H = gfd_prob{nvars:N0},
        ec_to_gecode_arith_expr1(E, H, 0, N0,N, [],Bs, Auxs0,AuxsT, GE, ConLev, Module),
        gfd_default_interval(Min, Max),
        update_vars_for_gecode(N0, N, Bs, H, Min, Max).


% convert existing gfd variable (or integer) to gecode representation
ec_to_gecode_oldvar(V, GV) :-
        ( integer(V) ->
            GV = V
        ; get_gecode_attr(V, gfd{idx:I,bool:B}) ->
            gfdvar(I, B, GV)
        ;
            fail
        ).


ec_to_gecode_var1(V, H, N0,N, GV) :-
        ( get_gecode_attr(V, Attr) ->
            % already converted to a gecode var 
            Attr = gfd{idx:I,bool:BI},
            N0 = N,
            gfdvar(I,BI, GV)
        ; 
            new_gfdvar(V, H, N0,N, GV) % convert to a new gecode var
        ).

ec_to_gecode_oldvar1(V, H, N0,N, OldGVs0,OldGVs, GV) :-
        ( get_gecode_attr(V, Attr) ->
            % already converted to a gecode var 
            Attr = gfd{idx:I,bool:BI},
            N0 = N,
            gfdvar(I,BI, GV),
            (I =< arg(nvars of gfd_prob, H) ->
                % already existing Gecode variable
                OldGVs = [GV|OldGVs0]
            ;
                OldGVs = OldGVs0
            )
        ; 
            OldGVs = OldGVs0,
            new_gfdvar(V, H, N0,N, GV) % convert to a new gecode var
        ).

ec_to_gecode_expr1(E, H, N0,N, Bs0,Bs, Auxs0,AuxsT, GE, ConLev, Module) :-
        boolean_expr(E), !,
        ec_to_gecode_bool_expr1(E, H, N0,N, Bs0,Bs, Auxs0,AuxsT, GE, ConLev, Module).
ec_to_gecode_expr1(E, H, N0,N, Bs0,Bs, Auxs0,AuxsT, GE, ConLev, Module) :-
        ec_to_gecode_arith_expr1(E, H, 0, N0,N, Bs0,Bs, Auxs0,AuxsT, GE, ConLev, Module).


ec_to_gecode_bool_expr1(V, H, N0,N, Bs0,Bs, Auxs0,AuxsT, GV, _ConLev, _Module) :- 
        var(V), !,
        Bs = [V|Bs0],
        Auxs0 = AuxsT,
        ec_to_gecode_var1(V, H, N0,N, GV).
ec_to_gecode_bool_expr1(subscript(T,S), H, N0,N, Bs0,Bs, Auxs0,AuxsT, GV, ConLev, Module) :- 
        subscript(T,S, E),
        ec_to_gecode_bool_expr1(E, H, N0,N, Bs0,Bs, Auxs0,AuxsT, GV, ConLev, Module).
ec_to_gecode_bool_expr1(1, _H, N0,N, Bs0,Bs, Auxs0,AuxsT, GE, _ConLev, _Module) :-
        N0 = N,
        Bs0 = Bs,
        Auxs0 = AuxsT,
        GE = 1.
ec_to_gecode_bool_expr1(0, _H, N0,N, Bs0,Bs, Auxs0,AuxsT, GE, _ConLev, _Module) :-
        N0 = N,
        Bs0 = Bs,
        Auxs0 = AuxsT,
        GE = 0.
ec_to_gecode_bool_expr1(eval(E), H, N0,N, Bs0,Bs, Auxs0,AuxsT, GV, ConLev, Module) :- 
        ec_to_gecode_bool_expr1(E, H, N0,N, Bs0,Bs, Auxs0,AuxsT, GV, ConLev, Module).
ec_to_gecode_bool_expr1(E, H, N0,N, Bs0,Bs, Auxs0,AuxsT, GE, ConLev, Module) :- 
        connective(E, Connective, SubExprs), !,
        ( foreach(SubE, SubExprs), 
          fromto(N0, N1,N2, N),
          fromto(Bs0, Bs1,Bs2, Bs),
          fromto(Auxs0, Auxs1,Auxs2, AuxsT),
          param(H,Module,ConLev),
          foreach(GSubE, GSubExprs)
        do
            ec_to_gecode_bool_expr1(SubE, H, N1,N2, Bs1,Bs2, Auxs1,Auxs2, GSubE, ConLev, Module)
        ),
        GE =.. [Connective|GSubExprs].
ec_to_gecode_bool_expr1(E, H, N0,N, Bs0,Bs, Auxs0,AuxsT, GE, ConLev, Module) :-
        % rel. constraints can be inlined in bool. expr,
        relation_constraint(E, Op, EX, EY), !, 
        ec_to_gecode_reifiedrc1(Op, EX, EY, H, N0,N, Bs0,Bs, Auxs0,AuxsT, GE, ConLev, Module).
ec_to_gecode_bool_expr1(C, H, N0,N, Bs0,Bs, Auxs0,AuxsT, GB, ConLev, Module) :-
        reifiable_constraint(C), !, 
        Auxs0 = [RCEvent|Auxs1],
        ec_to_gecode_reified1(C, H, N0,N, Bs0,Bs, Auxs1,AuxsT, RCEvent, GB, ConLev, Module).
ec_to_gecode_bool_expr1(C, _H, _N0,_N, _Bs0,_Bs, _Auxs0,_AuxsT, _GB, _ConLev, _Module) :-
        \+ integer(C),   % non-0/1 integer should fail, not raise error
        printf(error, "Boolean expression expected for: %w%n", [C]),
        set_bip_error(21).


ec_to_gecode_arith_expr1(V, H, _Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GV, _ConLev, _Module) :-
        var(V), !,
        Auxs = AuxsT,
        Bs = Bs0,
        ec_to_gecode_var1(V, H, N0,N, GV).
ec_to_gecode_arith_expr1(I, _H, _Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GI, _ConLev, _Module) :-
        integer(I), !,
        N0 = N,
        Bs0 = Bs,
        Auxs = AuxsT,
        GI = I.
ec_to_gecode_arith_expr1(eval(E), H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GE, ConLev, Module) ?- !,
        % ic compatibility: gfd expressions are always evaluated at runtime
        ec_to_gecode_arith_expr1(E, H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GE, ConLev, Module).
ec_to_gecode_arith_expr1(subscript(T,S), H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GE, ConLev, Module) ?- !,
        subscript(T,S,E),
        ec_to_gecode_arith_expr1(E, H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GE, ConLev, Module).
ec_to_gecode_arith_expr1(sum(Cs0*L0), H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GSP, ConLev, Module) ?- !,
        check_collection_to_list(flatten(L0),L), 
        check_collection_to_list(flatten(Cs0),Cs),
	(foreach(C, Cs) do check_integer(C)), 
        ec_to_gecode_arith_exprlist1(L, H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GL, ConLev, Module),
        GArr =.. [[]|GL],
	CArr =.. [[]|Cs],
	arity(GArr) =:= arity(CArr),
        GSP = sum(CArr,GArr).
ec_to_gecode_arith_expr1(sum(L0), H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GSum, ConLev, Module) ?- !,
        check_collection_to_list(flatten(L0),L),
	ec_to_gecode_arith_exprlist1(L, H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GL, ConLev, Module),
        GArr =.. [[]|GL],
        GSum = sum(GArr).
ec_to_gecode_arith_expr1(min(L0), H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GSum, ConLev, Module) ?- !,
        check_collection_to_list(flatten(L0),L), 
        ec_to_gecode_arith_exprlist1(L, H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GL, ConLev, Module),
        GArr =.. [[]|GL],
        GSum = min(GArr).
ec_to_gecode_arith_expr1(max(L0), H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GSum, ConLev, Module) ?- !,
        check_collection_to_list(flatten(L0),L), 
        ec_to_gecode_arith_exprlist1(L, H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GL, ConLev, Module),
        GArr =.. [[]|GL],
        GSum = max(GArr).
ec_to_gecode_arith_expr1(E0/E1, H, 0, N0,N, Bs0,Bs, Auxs,AuxsT, GV, ConLev, Module) ?- !,
        % IC compatible /
        new_gfdvar(_V, H, N0,N1, GV),
        ec_to_gecode_arith_expr1(E0, H, 0, N1,N2, Bs0,Bs1, Auxs1,Auxs2,
                                 GE0, ConLev, Module),
        ec_to_gecode_arith_expr1(E1, H, 0, N2,N, Bs1,Bs, Auxs2,AuxsT,
                                 GE1, ConLev, Module),
        Auxs = [post_rc(ConLev, (GE0 #= GE1*GV))|Auxs1]. 
ec_to_gecode_arith_expr1(sqrt(E0), H, 0, N0,N, Bs0,Bs, Auxs,AuxsT, GV, ConLev, Module) ?- !,
        % IC compatible sqrt
        new_gfdvar(_V, H, N0,N1, GV),
        ec_to_gecode_arith_expr1(E0, H, 0, N1,N, Bs0,Bs, Auxs1,AuxsT, GE0, ConLev, Module),
        Auxs = [post_rel(gfd_bc, GV, (#>=), 0), post_rc(ConLev, (GE0 #= GV*GV))|Auxs1]. 
ec_to_gecode_arith_expr1(element_g(Idx,Vs), H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GElm, ConLev, Module) ?- !,
        ec_to_gecode_arith_expr1(Idx, H, Lin, N0,N1, Bs0,Bs, Auxs,AuxsT, GIdx, ConLev, Module),
	check_collection_to_list(flatten(Vs), List),
        ec_to_gecode_varlist1(List, H, N1,N, GList, _),
	GArr =.. [[]|GList],
        % order of arg in Gecode is reverse of ECLiPSe 
	GElm = element(GArr, GIdx).
ec_to_gecode_arith_expr1(E0, H, Lin, N0,N, Bs0,Bs, Auxs,AuxsT, GE, ConLev, Module) :-
        inline_op(E0, E), !,
        functor(E, Name, Arity),
        functor(GE, Name, Arity),
        ( foreacharg(Arg, E), foreacharg(GArg, GE),
          fromto(N0, N1,N2, N), 
          fromto(Bs0, Bs1,Bs2, Bs),
          param(H,Module,Lin,ConLev),
          fromto(Auxs, Auxs1,Auxs2, AuxsT)
        do
            ec_to_gecode_arith_expr1(Arg, H, Lin, N1,N2, Bs1,Bs2, Auxs1,Auxs2,  GArg, ConLev, Module)
        ).
% following clauses are for expressions that needs to be factored out as aux. 
% these used to be the non-linear parts but Gecode can now handle most
% non-linear parts directly
ec_to_gecode_arith_expr1(E, H, 0, N0,N, Bs0,Bs, Auxs0,AuxsT, GNewV, ConLev, Module) :-
        aux_op(E, EventTemp, NewV, GNewV, ArgType, ConLev), !,
        new_gfdvar(NewV, H, N0,N1, GNewV),
        ec_to_gecode_aux_op1(ArgType, E, H, N1,N, Bs0,Bs, Auxs0,AuxsT, EventTemp, ConLev, Module).
ec_to_gecode_arith_expr1(C, H, 0, N0,N, Bs0,Bs, Auxs0,AuxsT, GB, ConLev, Module) :-
        reifiable_constraint(C), !, 
        Auxs0 = [RCEvent|Auxs1],
        ec_to_gecode_reified1(C, H, N0,N, Bs0,Bs, Auxs1,AuxsT, RCEvent, GB, ConLev, Module).
ec_to_gecode_arith_expr1(E, _H, In, _N0,_N, _Bs0,_Bs, _Auxs0,_AuxsT,
                         _GNewV, _ConLev, _Module) :-
        write(error, "Unrecognised constraint subexpression "),
        (In == 1 -> 
            write(error, "(non-inlined subexpressions not allowed) ")
        ;
            true
        ),
        printf(error, "in: %w%n", [E]),
        set_bip_error(21).


%------------------------------------------------------------------------
% Constraints


'::_body'(X, Domain, Module):-
        get_prob_handle_nvars(H, NV0),
        normalise_vars(X, NX),
	process_domain_domain(Domain, NormalDomain, Module),
        process_domain_vars(NX, NormalDomain, H, NV0,NV, [],OldGVs, _),
        !, % cut here to avoid posting events with a live choicepoint
        assign_domain(NormalDomain, H, NV, OldGVs).
'::_body'(X, Domain, _Module) :-
        get_bip_error(E),
        error(E,(X :: Domain)).


:- tool('::'/3, '::_body'/4).
:- tool('#::'/3, '::_body'/4).


'::_body'(X, Domain, Bool, Module):-
        get_prob_handle_nvars(H, N0),
        gfd_default_interval(Min, Max),
        ec_to_gecode_domain_reified1(X, Domain, Bool, H, N0,N, [],Bs, Event, _, Module),
        !,
        update_vars_for_gecode(N0, N, Bs, H, Min, Max),
        post_new_event(Event, H).
'::_body'(X, Domain, Bool, _Module):-
        get_bip_error(E),
        error(E, ::(X, Domain, Bool)).


ec_to_gecode_domain_reified1(X, Domain, Bool, H, N0,N, Bs0,Bs, Event, GBool, Module) :-
        process_domain_domain(Domain, NDomain, Module),
        ( var(Bool) ->
            ec_to_gecode_var1(Bool, H, N0,N1, GBool),
            Bs = [Bool|Bs0]
        ; integer(Bool) ->
            Bool >= 0, 
            Bool =< 1, % fails if Bool out of range (ic compatible)
            GBool = Bool,
            N0 = N1,
            Bs = Bs0
        ;
            set_bip_error(5)
        ),
        ( var(X) ->
            ec_to_gecode_var1(X, H, N1, N, GX)
        ; integer(X) ->
            GX = X,
            N1 = N
        ; 

            set_bip_error(5)
        ),
        % domain is normalised and type checked already
        ( NDomain = [I..I] ->
            Event = post_var_val_reif(GX, I, GBool)
        ; 
            NDomain = [Lo..Hi|T],
            ( T == [] ->
                % Domain is a simple interval
                Event = post_var_interval_reif(GX,Lo,Hi, GBool)
            ;
                DArray =.. [[]|NDomain],
                Event = post_var_dom_reif(GX,DArray,GBool)
            )
        ).


normalise_vars(V, N) :-
        var(V), !,
        N = [V].
normalise_vars(I, N) :-
        integer(I), !,
        N = [I].
normalise_vars(Xs, NXs) :-
        (check_collection_to_list(flatten(Xs), NXs) -> true ; set_bip_error(5)).


split_first_domain([H0|T0], H, T) ?- !,
        H0 = H, T0 = T.
split_first_domain(Dom, H, T) :-
        ( nonvar(Dom) ->
            H = Dom,
            T = []
        ;
            set_bip_error(4)
        ).

process_domain_domain(Domain, NormalDomain, Module) :- 
        split_first_domain(Domain, H, T),
        subdomain(H, Lo, Hi, Module),
	( T \== [] ->
	    ( Lo =< Hi ->
		Domain1 = [Lo..Hi | Domain0]
	    ;
		Domain1 = Domain0
	    ),
	    (
		foreach(Sub, T),
		fromto(Domain0, Out, In, []),
		param(Module)
	    do
		subdomain(Sub, Lo, Hi, Module),
		% Filter empty ranges (Lo > Hi).
		( Lo =< Hi ->
		    Out = [Lo..Hi | In]
		;
		    Out = In
		)
	    ),
	    % Order the intervals.
	    number_sort(2, =<, Domain1, SortedUpperBoundsDomain),
	    number_sort(1, =<, SortedUpperBoundsDomain, SortedIntervalDomain),
	    [Lo0..Hi0 | SortedRest] = SortedIntervalDomain,
	    % Collapse zero width intervals to constants and merge
	    % overlapping/adjacent subdomains.  
	    (
		foreach(Lo..Hi, SortedRest),
		fromto(Lo0..Hi0, LoIn..HiIn, LoOut..HiOut, FinalSubDomain),
		fromto(NormalDomain, In, Out, [FinalSubDomain])
	    do
		( HiIn + 1 >= Lo ->
		    % There is no gap between HiIn and Lo so merge
		    In = Out,
		    LoOut = LoIn,
		    HiOut is max(Hi, HiIn)
		;
		    % There is a gap between HiIn and Lo
		    In = [LoIn..HiIn | Out],
		    LoOut = Lo,
		    HiOut = Hi
		)
	    )
	;
	    NormalDomain = [Lo..Hi]
	).

bound(I, B, _Module) :-
        integer(I), !,
        B = I.
bound(A, B, Module) :-
%        ground(A),
        % cut after subcall in case A is non-determinate
        block(subcall((B is A, !), [])@Module, _Tag, set_bip_error(21)),
        check_integer(B).


subdomain(Lo..Hi, Lo1, Hi1, Module) ?- !,
        bound(Lo, Lo1, Module),
        bound(Hi, Hi1, Module).
subdomain(I, Lo1, Hi1, Module) :-
        Hi1 = Lo1,
        bound(I, Lo1, Module).

process_domain_vars([V1|Vs], Domain, H, NV0,NV, OldGVs0,OldGVs, [GV1|GVs]) :-
        var(V1), !,
        ( Domain = [B..B], \+ ec_to_gecode_oldvar(V1,_) ->  
            % singleton domain, and a nnon-gecode variable
            V1 = B, % just assign it
            GV1 = B,
            NV0 = NV1,
            OldGVs0 = OldGVs1
        ;
            ec_to_gecode_oldvar1(V1, H, NV0,NV1, OldGVs0,OldGVs1, GV1)
        ),
        process_domain_vars(Vs, Domain, H, NV1,NV, OldGVs1,OldGVs, GVs).
process_domain_vars([I|Vs], Domain, H, NV0,NV, OldGVs0,OldGVs, [I|GVs]) :-
        integer(I), !,
        is_in_given_domain(I, Domain),
        process_domain_vars(Vs, Domain, H, NV0,NV, OldGVs0,OldGVs, GVs).
process_domain_vars([_|_], _,_,_,_,_,_,_) :- !,
        set_bip_error(5).
process_domain_vars([], _D, _H, NV,NV, OldGVs,OldGVs, []).


is_in_given_domain(I, [Lo..Hi|Ds]) :-
        (I >= Lo, I =< Hi -> true ; is_in_given_domain(I, Ds)).


assign_domain(Domain, H, NV, OldGVs) :- 
        Domain = [Lo..Hi|T], !,
        (T == [] ->
            Hi >= Lo,
            assign_domain_interval(H, NV, OldGVs, Lo, Hi)
        ;
            assign_multi_domain_intervals(H, NV, OldGVs, Domain)
        ).

assign_domain_interval(H, NV, OldGVs, Lo, Hi) :-
        ( OldGVs == [] ->
            true
        ;
            GVArr =.. [[]|OldGVs],
            post_new_event(post_interval(GVArr,Lo,Hi), H)
        ),
        update_newvars_with_domain_interval(H, NV, Lo, Hi).

assign_domain_interval1(H, NV0, NV, OldGVs, Lo, Hi) :-
        ( OldGVs == [] ->
            true
        ;
            GVArr =.. [[]|OldGVs],
            post_new_event(post_interval(GVArr,Lo,Hi), H)
        ),
        ( NV0 == NV ->
            true
        ;
            do_update_newvars_with_domain_interval(H, NV, Lo, Hi)
        ).

update_newvars_with_domain_interval(H, NV, Lo, Hi) :-
        ( NV =:= arg(nvars of gfd_prob, H) -> % no new vars 
            true
        ;
            do_update_newvars_with_domain_interval(H, NV, Lo, Hi)
        ).

do_update_newvars_with_domain_interval(H, NV, Lo, Hi) :-
        setarg(nvars of gfd_prob, H, NV),
        post_new_event(newvars_interval(NV,Lo,Hi), H).

assign_multi_domain_intervals(H, NV, OldGVs, Domain) :-
        DArray =.. [[]|Domain],
        ( OldGVs == [] ->
            true
        ;
            GVArr =.. [[]|OldGVs],
            post_new_event(post_dom(GVArr,DArray), H)
        ),
        update_newvars_with_multi_domain_intervals(H, NV, DArray).

update_newvars_with_multi_domain_intervals(H, NV, DArray) :-
        ( NV =:= arg(nvars of gfd_prob, H) -> % no new vars 
            true
        ;
            setarg(nvars of gfd_prob, H, NV),
            post_new_event(newvars_dom(NV,DArray), H)
        ).



alldifferent(Vars) :-
        alldifferent_c(Vars, default).


alldifferent_c(Vars, ConLev) :-
        check_collection_to_list(Vars, List),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1(List, H, NV0,NV, GList, _),
        GArray =.. [[]|GList], !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_alldiff(ConLev,GArray), H).
alldifferent_c(Vars, _ConLev) :-
        get_bip_error(E),
        error(E, alldifferent(Vars)).


alldifferent_cst(Vars, Offsets) :-
        alldifferent_cst_c(Vars, Offsets, default).


alldifferent_cst_c(Vars, Offsets, ConLev) :-
        check_collection_to_list(Vars, List),
        check_collection_to_list(Offsets, OffList),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1(List, H, NV0,NV, GList, _),
        (foreach(O,  OffList) do check_integer(O)),
        GArray =.. [[]|GList],
        OArray =.. [[]|OffList],
        arity(GArray, N),
        (arity(OArray, N) -> true ; set_bip_error(6)), !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_alldiff_offsets(ConLev,GArray,OArray), H).
alldifferent_cst_c(Vars, Offsets, _ConLev) :-
        get_bip_error(E),
        error(E, alldifferent_cst(Vars,Offsets)).


nvalues(Vars, Rel0, N) :-
        is_valid_rel_op(Rel0, Rel),
        check_collection_to_list(Vars, VList),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([N|VList], H, NV0,NV, [GN|GVars], _),
        VArr =.. [[]|GVars], !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_nvalues(VArr, Rel, GN), H).
nvalues(Vars, Rel0, N) :-
        get_bip_error(E),
        error(E, nvalues(Vars, Rel0, N)).


mem(Vars, X) :-
        mem_c(Vars, X, default).

mem_c(Vars, X, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        mem_events1(Vars, X, ConLev, H, NV0,NV, Event,EventsT),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event_with_aux(Event,EventsT, H).
mem_c(Vars, X, _ConLev) :-
        get_bip_error(E),
        error(E, mem(Vars, X)).

mem_events1(Vars, X, _ConLev, H, NV0,NV, Event,EventsT) :-
        check_collection_to_list(Vars, VList),
        ec_to_gecode_varlist1([X|VList], H, NV0,NV, [GX|GVars], _),
        VArr =.. [[]|GVars],
        Event = [post_mem(VArr, GX)|EventsT].


mem(Vars, X, Bool) :-
        mem_reif_c(Vars, X, Bool, default).

mem_reif_c(Vars, X, Bool, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        mem_reif_event1(Vars, X, Bool, ConLev, H, NV0,NV, [],Bs, Event, _GBool),
        !,
        gfd_default_interval(Min, Max),
	update_vars_for_gecode(NV0, NV, Bs, H, Min, Max),
        post_new_event_with_aux([Event|EventsT], EventsT, H).
mem_reif_c(Vars, X, Bool, _ConLev) :-
        get_bip_error(E),
        error(E, mem(Vars, X, Bool)).
        
mem_reif_event1(Vars, X, Bool, _ConLev, H, NV0,NV, Bs0,Bs, Event, GBool) :-
        check_collection_to_list(Vars, VList),
        ec_to_gecode_varlist1([X|VList], H, NV0,NV1, [GX|GVars], _),
        ( var(Bool) ->
            ec_to_gecode_var1(Bool, H, NV1,NV, GBool),
            Bs = [Bool|Bs0]
        ; integer(Bool) ->
            Bool >= 0,
            Bool =< 1,
            GBool = Bool,
            NV1 = NV,
            Bs = Bs0
        ;
            set_bip_error(5)
        ),
        VArr =.. [[]|GVars],
        Event = post_mem_reif(VArr, GX, GBool).


count(Value, Vars, Rel, N) :-
        count_c(Value, Vars, Rel, N, default).

count_c(Value, Vars, Rel, N, ConLev) :-
	get_prob_handle_nvars(H, NV0),
	count_events1(Value, Vars, Rel, N, ConLev, H, NV0,NV, CountEvents, EventsT),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event_with_aux(CountEvents, EventsT, H).
count_c(Value, Vars, Rel, N, _ConLev) :-
        get_bip_error(E),
        error(E, count(Value, Vars, Rel, N)).

count_events1(Value, Vars, Rel0, N, ConLev, H, NV0,NV, Events, EventsT) :-
        is_valid_rel_op(Rel0, Rel),
        check_collection_to_list(Vars, List),
        ec_to_gecode_varlist1([Value,N|List], H, NV0,NV, [GValue,GN|GList],_),
        GArray =.. [[]|GList],
	CEvent = post_count(ConLev,GValue,GArray,Rel,GN), 
        ( integer(N) ->
            N =< arity(GArray),
	    Events = [CEvent|EventsT]
        ; var(N) ->
	    Hi is arity(GArray),
	    % GN will be in Gecode by the time this event is executed
	    Events = [post_interval([](GN), 0, Hi),CEvent|EventsT] 
        ;
            set_bip_error(5)
        ).


% compatibility
occurrences(Value, Vars, N) :-
        count_c(Value, Vars, '#=', N, default).

occurrences_c(Value, Vars, N, ConLev) :-
        count_c(Value, Vars, '#=', N, ConLev).

% compatibility
atmost(N, Vars, Value) :-
        count_c(Value, Vars, '#=<', N, default).

atmost_c(N, Vars, Value, ConLev) :-
        count_c(Value, Vars, '#=<', N, ConLev).

atleast(N, Vars, Value) :-
        count_c(Value, Vars, '#>=', N, default).

atleast_c(N, Vars, Value, ConLev) :-
        count_c(Value, Vars, '#>=', N, ConLev).


among_body(Values, Vars, Rel, N, Module) :-
        among_c(Values, Vars, Rel, N, default, Module).

among_c(Values, Vars, Rel, N, ConLev, Module) :-
        get_prob_handle_nvars(H, NV0),
        among_events1(Values, Vars, Rel, N, ConLev, H, NV0,NV, 
                     Event,EventsT, Module),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event_with_aux(Event,EventsT, H).
among_c(Values, Vars, Rel, N, _ConLev, _Module) :-
        get_bip_error(E),
        error(E, among(Values, Vars, Rel, N)).

among_events1(ValSpec, Vars, Rel0, N, _ConLev, H, NV0,NV, Ev,EvT, Module) :-
        is_valid_rel_op(Rel0, Rel),
        check_collection_to_list(ValSpec, SpecList),
        ( foreach(Spec,SpecList), param(Module),
          foreach(Lo..Hi, Vals)
        do
            subdomain(Spec, Lo, Hi, Module)
        ),
        ValsArr =.. [[]|Vals],
        check_collection_to_list(Vars, VarsL),
        ec_to_gecode_varlist1([N|VarsL], H, NV0, NV, [GN|GVars], _),
        VarsArr =.. [[]|GVars],
        % gac only currently, no need for ConLev
        Ev = [post_among(ValsArr, VarsArr, Rel, GN)|EvT].
                       

count_matches(Values, Vars, Rel, N) :-
        count_matches_c(Values, Vars, Rel, N, default).

count_matches_c(Values, Vars, Rel, N, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        count_matches_events1(Values, Vars, Rel, N, ConLev, H, NV0,NV, 
                     Event,EventsT),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event_with_aux(Event,EventsT, H).
count_matches_c(Values, Vars, Rel, N, _ConLev) :-
        get_bip_error(E),
        error(E, count_matches(Values, Vars, Rel, N)).

count_matches_events1(Vals, Vars, Rel0, N, _ConLev, H, NV0,NV, Ev,EvT) :-
        is_valid_rel_op(Rel0, Rel),
        check_collection_to_list(Vals, ValsL),
        ValsArr =.. [[]|ValsL],
        check_collection_to_list(Vars, VarsL),
        ec_to_gecode_varlist1([N|VarsL], H, NV0, NV, [GN|GVars], _),
        VarsArr =.. [[]|GVars],
        (arity(VarsArr) =:= arity(ValsArr) -> true ; set_bip_error(6)),
        % gac only currently, no need for ConLev
        Ev = [post_count_matches(ValsArr, VarsArr, Rel, GN)|EvT].
                       

element(Index, Collection, Value) :-
        element_body(Index, Collection, Value, ecl, default).

element_c(Index, Collection, Value, ConLev) :-
        element_body(Index, Collection, Value, ecl, ConLev).

element_g(Index, Collection, Value) :-
        element_body(Index, Collection, Value, gc, default).

element_g_c(Index, Collection, Value, ConLev) :-
        element_body(Index, Collection, Value, gc, ConLev).

element_body(Index, Collection, Value, IndexType, ConLev) :-
	get_prob_handle_nvars(H, NV0),
	element_events1(Index, Collection, Value, IndexType, ConLev, H, NV0,NV, Events, EventsT),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event_with_aux(Events, EventsT, H).
element_body(Index, Collection, Value, IndexType, _ConLev) :-
        get_bip_error(E),
        ( IndexType == ecl ->
            error(E,element(Index, Collection, Value))
        ;
            error(E,element_g(Index, Collection, Value))
        ).

element_events1(Index, Collection, Value, IndexType, ConLev, H, NV0,NV, Es, EsT) :-
        check_collection_to_list(Collection, List),
        ec_to_gecode_varlist1([Index,Value|List], H, NV0,NV, [GIndex,GValue|GList],_),
	( IndexType == ecl ->
	    Array =.. [[],0|GList],  % add a dummy first element for index 0
	    Lo = 1
	;
	    Array =.. [[]|GList],
	    Lo = 0
        ),
        Hi is arity(Array)-1,
        EEvent = post_element(ConLev, GIndex, Array, GValue),  
	( integer(Index) ->
	   Index >= Lo,
	   Index =< Hi,
	   Es = [EEvent|EsT] 
	;
	   
	   Es = [post_interval([](GIndex), Lo, Hi),EEvent|EsT]
	).

:- export struct(gcc(low,high,value)),
          struct(occ(occ,value)).


gcc(BoundsList, Vars) :-
        gcc_c(BoundsList, Vars, default).

gcc_c(BoundsList, Vars, ConLev) :-        
        check_collection_to_list(Vars, VList),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1(VList, H, NV0,NV1, GVList, _),
        GVs =.. [[]|GVList],
        Bounds =.. [[]|BoundsList],
        arity(Bounds, M),
        dim(Vals, [M]),
        dim(Occurrences, [M]), 
        ( foreacharg(Spec, Bounds),
          fromto(NV1, NV2,NV3, NV),
          fromto(Ev, Ev1,Ev2, EvT0),
          param(H),
          foreacharg(Value, Vals),
          foreacharg(GOcc, Occurrences)
        do
            translate_gcc_spec_events1(Spec, Value, H, NV2,NV3, Ev1,Ev2, GOcc)
        ), !,
        update_gecode_with_default_newvars(H, NV0, NV1),
        setarg(nvars of gfd_prob, H, NV),
        % the gcc event must be the last event posted here
        EvT0 = [post_gcc(ConLev, Vals, Occurrences, GVs)|EvT],
        post_new_event_with_aux(Ev, EvT, H).
gcc_c(BoundsList, Vars, _ConLev) :-
        get_bip_error(Error),
        error(Error, gcc(BoundsList, Vars)).

translate_gcc_spec_events1(gcc{low:Lo,high:Hi,value:Val0}, Val, H, NV0,NV,
                           Ev0,EvT, GOcc) ?- !,
        check_integer(Val0),
        Val0 = Val,
        check_integer(Hi),
        check_integer(Lo),
        new_gfdvar(_Occ, H, NV0,NV, GOcc),
        Ev0 = [newvars_interval(NV, Lo, Hi)|EvT].
translate_gcc_spec_events1(occ{occ:Occ0,value:Val0}, Val, _H, NV0,NV,
                           Ev0,EvT, GOcc) ?- !,
        check_integer(Val0),
        Val0 = Val,
        NV0 = NV,
        Ev0 = EvT,
        ( ec_to_gecode_oldvar(Occ0, GOcc) ->
            true
        ;
            set_bip_error(5) % not a gfd domain var or integer
        ).
translate_gcc_spec_events1(_Spec, _Val, _H, _,_, _,_, _Occ) :-
        set_bip_error(6).  % range error -- spec not recognised


sorted(Us0, Ss0) :-
        sorted_c(Us0, Ss0, default).

sorted_c(Us0, Ss0, ConLev) :-
        ( var(Us0) -> 
            Us0 = Us
        ;
            check_collection_to_list(Us0, Us)
        ),
        ( var(Ss0) -> 
            nonvar(Us0),
            Ss0 = Ss
        ;
            check_collection_to_list(Ss0, Ss)
        ),
	( foreach(_,Us), foreach(_,Ss) do true ),
        !,
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1(Ss, H, NV0,NV1, GSs, _),
        SsArr =.. [[]|GSs],
        ec_to_gecode_varlist1(Us, H, NV1,NV, GUs, _),
        UsArr =.. [[]|GUs], 
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_sorted2(ConLev, UsArr, SsArr), H).
sorted_c(Us0, Ss0, _ConLev) :-
        get_bip_error(E),
        error(E, sorted(Us0, Ss0)).


sorted(Us0, Ss0, Ps0) :-
        sorted_body(Us0, Ss0, Ps0, ecl, default).

sorted_c(Us0, Ss0, Ps0, ConLev) :-
        sorted_body(Us0, Ss0, Ps0, ecl, ConLev).

sorted_g(Us0, Ss0, Ps0) :-
        sorted_body(Us0, Ss0, Ps0, gc, default).

sorted_g_c(Us0, Ss0, Ps0, ConLev) :-
        sorted_body(Us0, Ss0, Ps0, gc, ConLev).

sorted_body(Us0, Ss0, Ps0, IndexType, ConLev) :-
        ( var(Us0) -> 
            Us0 = Us
        ;
            check_collection_to_list(Us0, Us)
        ),
        ( var(Ss0) -> 
            Ss0 = Ss
        ;
            check_collection_to_list(Ss0, Ss)
        ),
        ( var(Ps0) -> 
            Ps0 = Ps
        ;
            check_collection_to_list(Ps0, Ps)
        ),
	get_prob_handle_nvars(H, NV0),
        convert_sorted_lists1(IndexType, Us, Ss, Ps, H, NV0,NV1,NV, PMin, PMax,  
                              UsArr, SsArr, PsArr, OldPs), 
        !,
        update_gecode_with_default_newvars(H, NV0, NV1),
        assign_domain([PMin..PMax], H, NV, OldPs),
        post_new_event(post_sorted(ConLev, UsArr, SsArr, PsArr), H).
sorted_body(Us0, Ss0, Ps0, IndexType, _ConLev) :-
        get_bip_error(E),
        (IndexType = ecl ->
            error(E, sorted(Us0, Ss0, Ps0))
        ;
            error(E, sorted_g(Us0, Ss0, Ps0))
        ).

convert_sorted_lists1(ecl, Us, Ss, Ps, H, NV0,NV1,NV, 1, PMax,  UsArr, SsArr, PsArr, OldPs) :-
        gfd_maxint(Max),
	( foreach(U,Us), foreach(_,Ss), foreach(_,Ps),
          fromto(Max, Min0,Min1, Min), count(_, 1, PMax),
          foreach(GU, GUs)
        do 
           get_min(U, UMin),
           (UMin < Min0 -> Min1 is UMin - 1 ; Min1 = Min0),
           ( var(U) ->
               % must be a gfd var because of get_min/2
               get_gecode_var(U, GU)
           ; integer(U) ->
               GU = U
           ;
               set_bip_error(5)
           )
        ),
        ec_to_gecode_varlist1(Ss, H, NV0,NV1, GSs, _),
        SsArr =.. [[],Min|GSs],
        UsArr =.. [[],Min|GUs],
        process_domain_vars(Ps, [1..PMax], H, NV1,NV, [],OldPs, GPs),
        PsArr =.. [[],0|GPs].
 convert_sorted_lists1(gc, Us, Ss, Ps, H, NV0,NV1,NV, 0, Max, UsArr, SsArr, PsArr, OldPs) :-
	( foreach(_,Us), foreach(_,Ss), foreach(_,Ps), 
          count(_,0, Max)
        do 
            true
        ),
        ec_to_gecode_multivarlists1([Ss,Us], H, NV0,NV1, [GSs,GUs]),
        SsArr =.. [[]|GSs],
        UsArr =.. [[]|GUs],
        process_domain_vars(Ps, [0..Max], H, NV1,NV, [],OldPs, GPs),
        PsArr =.. [[]|GPs].


circuit(Succ) :-
        circuit_offset_body(Succ, 0, ecl, default).

circuit_g(Succ) :-
        circuit_offset_body(Succ, 0, gc, default).

circuit_c(Succ, ConLev) :-
        circuit_offset_body(Succ, 0, ecl, ConLev).

circuit_g_c(Succ, ConLev) :-
        circuit_offset_body(Succ, 0, gc, ConLev).


circuit_offset(Succ, Off) :-
        circuit_offset_body(Succ, Off, ecl, default).

circuit_offset_g(Succ, Off) :-
        circuit_offset_body(Succ, Off, gc, default).

circuit_offset_c(Succ, Off, ConLev) :-
        circuit_offset_body(Succ, Off, ecl, ConLev).

circuit_offset_g_c(Succ, Off, ConLev) :-
        circuit_offset_body(Succ, Off, gc, ConLev).

circuit_offset_body(Succ, Off0, IndexType, ConLev) :-
        check_integer(Off0),
        check_collection_to_list(Succ, SList),
        length(SList, N), 
        (IndexType == ecl ->
            Off is Off0 + 1, Min = Off, Max is N + Off0 
        ;
            Off is Off0, Min = Off, Max is N + Off - 1
        ),
        get_prob_handle_nvars(H, NV0),
        process_domain_vars(SList, [Min..Max], H, NV0,NV, [],OldVs, GSs),
        !,
        assign_domain([Min..Max], H, NV, OldVs),
        SArr =.. [[]|GSs],
        post_new_event(post_circuit(ConLev,SArr,Off), H).
circuit_offset_body(Succ, Off0, IndexType, _ConLev) :-
        get_bip_error(E),
        (IndexType = ecl ->
            error(E, circuit_offset(Succ, Off0))
        ;
            error(E, circuit_offset_g(Succ, Off0))
        ).


circuit(Succ, CostMatrix, Cost) :-
        circuit_offset_body(Succ, 0, CostMatrix, [], Cost, ecl, default).

circuit(Succ, CostMatrix, ArcCosts, Cost) :-
        circuit_offset_body(Succ, 0, CostMatrix, ArcCosts, Cost, ecl, default).

circuit_g(Succ, CostMatrix, Cost) :-
        circuit_offset_body(Succ, 0, CostMatrix, [], Cost, gc, default).

circuit_g(Succ, CostMatrix, ArcCosts, Cost) :-
        circuit_offset_body(Succ, 0, CostMatrix, ArcCosts, Cost, gc, default).

circuit_c(Succ, CostMatrix, Cost, ConLev) :-
        circuit_offset_body(Succ, 0, CostMatrix, [], Cost, ecl, ConLev).

circuit_c(Succ, CostMatrix, ArcCosts, Cost, ConLev) :-
        circuit_offset_body(Succ, 0, CostMatrix, ArcCosts, Cost, ecl, ConLev).

circuit_g_c(Succ, CostMatrix, Cost, ConLev) :-
        circuit_offset_body(Succ, 0, CostMatrix, [], Cost, gc, ConLev).

circuit_g_c(Succ, CostMatrix, ArcCosts, Cost, ConLev) :-
        circuit_offset_body(Succ, 0, CostMatrix, ArcCosts, Cost, gc, ConLev).


circuit_offset(Succ, Off, CostMatrix, Cost) :-
        circuit_offset_body(Succ, Off, CostMatrix, [], Cost, ecl, default).

circuit_offset(Succ, Off, CostMatrix, ArcCosts, Cost) :-
        circuit_offset_body(Succ, Off, CostMatrix, ArcCosts, Cost, ecl, default).

circuit_offset_g(Succ, Off, CostMatrix, Cost) :-
        circuit_offset_body(Succ, Off, CostMatrix, [], Cost, gc, default).

circuit_offset_g(Succ, Off, CostMatrix, ArcCosts, Cost) :-
        circuit_offset_body(Succ, Off, CostMatrix, ArcCosts, Cost, gc, default).

circuit_offset_c(Succ, Off, CostMatrix, Cost, ConLev) :-
        circuit_offset_body(Succ, Off, CostMatrix, [], Cost, ecl, ConLev).

circuit_offset_c(Succ, Off, CostMatrix, ArcCosts, Cost, ConLev) :-
        circuit_offset_body(Succ, Off, CostMatrix, ArcCosts, Cost, ecl, ConLev).

circuit_offset_g_c(Succ, Off, CostMatrix, Cost, ConLev) :-
        circuit_offset_body(Succ, Off, CostMatrix, [], Cost, gc, ConLev).

circuit_offset_g_c(Succ, Off, CostMatrix, ArcCosts, Cost, ConLev) :-
        circuit_offset_body(Succ, Off, CostMatrix, ArcCosts, Cost, gc, ConLev).

circuit_offset_body(Succ, Off0, CostMatrix, ArcCosts, Cost, IndexType, ConLev) :-
        check_integer(Off0),
        check_collection_to_list(Succ, SList),
        length(SList, N), 
        (IndexType == ecl ->
            Off is Off0 + 1, Min = Off, Max is N + Off0 
        ;
            Off is Off0, Min = Off, Max is N + Off - 1
        ),
        get_prob_handle_nvars(H, NV0),
        process_domain_vars(SList, [Min..Max], H, NV0,NV1, [],OldVs, GSs),
        SArr =.. [[]|GSs],
        circuit_offset_events1(SArr, Off, CostMatrix, ArcCosts, Cost,
                               ConLev, H, NV1,NV, Event, EventTail),
        !,
        assign_domain([Min..Max], H, NV1, OldVs),
        % NV1 because assign_domain updated nvars to NV1
        update_gecode_with_default_newvars(H, NV1, NV),
        post_new_event_with_aux(Event, EventTail, H).
circuit_offset_body(Succ, Off0, CostMatrix, ArcCosts, Cost, IndexType, _ConLev) :-
        get_bip_error(E),
        (IndexType = ecl ->
            ( ArcCosts == [] ->
                error(E, circuit_offset(Succ, Off0, CostMatrix, Cost))
            ;
                error(E, circuit_offset(Succ, Off0, CostMatrix, ArcCosts, Cost))
            )
        ;
            ( ArcCosts == [] ->
                error(E, circuit_offset_g(Succ, Off0, CostMatrix, Cost))
            ;
                error(E, circuit_offset_g(Succ, Off0, CostMatrix, ArcCosts, Cost))
            )
        ).

circuit_offset_events1(SArr, Off, CostMatrix, ArcCosts, Cost, ConLev, H,
                       NV1,NV, Event, EventTail) :-
        check_collection_to_list(flatten(CostMatrix), CMList),
        CMArr =.. [[]|CMList],
        ec_to_gecode_var1(Cost, H, NV1,NV2, GCost),
        ( ArcCosts == [] ->
            NV = NV2,
            GAC = []
        ;
            check_collection_to_list(ArcCosts, ACList),
            ec_to_gecode_varlist1(ACList, H, NV2,NV, GACList, _),
            GAC =.. [[]|GACList]
        ),
        Event = [post_circuit_cost(ConLev,SArr,CMArr,GAC,GCost,Off)|EventTail].

% The common code for posting a aux. circuit (with cost) constraint in a
% constraint expression, we cannot treat the Succ variables differently here
circuit_aux_events1(Succ, Off, CostMat, ArcCosts, Cost, ConLev, H, N0,N, Auxs0,AuxsT) :-
        check_collection_to_list(Succ, SuccL),
        ec_to_gecode_varlist1(SuccL, H, N0,N1, GSuccL, _),
        SArr =.. [[]|GSuccL],
        circuit_offset_events1(SArr, Off, CostMat, ArcCosts, Cost, ConLev, H, N1,N, Auxs0,AuxsT). 


ham_path(Start, End, Succ) :-
        ham_path_offset_body(Start, End, Succ, 0, ecl, default).

ham_path_g(Start, End, Succ) :-
        ham_path_offset_body(Start, End, Succ, 0, gc, default).

ham_path_c(Start, End, Succ, ConLev) :-
        ham_path_offset_body(Start, End, Succ, 0, ecl, ConLev).

ham_path_g_c(Start, End, Succ, ConLev) :-
        ham_path_offset_body(Start, End, Succ, 0, gc, ConLev).


ham_path_offset(Start, End, Succ, Off) :-
        ham_path_offset_body(Start, End, Succ, Off, ecl, default).

ham_path_offset_g(Start, End, Succ, Off) :-
        ham_path_offset_body(Start, End, Succ, Off, gc, default).

ham_path_offset_c(Start, End, Succ, Off, ConLev) :-
        ham_path_offset_body(Start, End, Succ, Off, ecl, ConLev).

ham_path_offset_g_c(Start, End, Succ, Off, ConLev) :-
        ham_path_offset_body(Start, End, Succ, Off, gc, ConLev).


ham_path_offset_body(Start, End, Succ, Off0, IdxType, ConLev) :-
        check_integer(Off0),
        check_collection_to_list(Succ, SList),
        length(SList, Size),
        % The successor to End node is a dummy id = largest id + 1 (Max)
        ( IdxType == ecl ->
            Off is Off0 + 1, Min = Off, Max is Size + Off
        ;
            Off is Off0, Min = Off, Max is Size + Off
        ),
        get_prob_handle_nvars(H, NV0),
        process_domain_vars([Start,End|SList], [Min..Max], H, NV0,NV, 
                            [],OldVs, [GStart,GEnd|GSs]),
        !,
        assign_domain([Min..Max], H, NV, OldVs),
        SArr =.. [[]|GSs],
        post_new_event(post_ham_path(ConLev,GStart,GEnd,SArr,Off), H).
ham_path_offset_body(Start, End, Succ, Off0, IdxType, _ConLev) :-
        get_bip_error(E),
        (IdxType = ecl ->
            error(E, ham_path_offset(Start, End, Succ, Off0))
        ;
            error(E, ham_path_offset(Start, End, Succ, Off0))
        ).


ham_path(Start, End, Succ, CostMatrix, Cost) :-
        ham_path_offset_body(Start, End, Succ, 0, CostMatrix, [], Cost, ecl, default).

ham_path_g(Start, End, Succ, CostMatrix, Cost) :-
        ham_path_offset_body(Start, End, Succ, 0, CostMatrix, [], Cost, gc, default).

ham_path_c(Start, End, Succ, CostMatrix, Cost, ConLev) :-
        ham_path_offset_body(Start, End, Succ, 0, CostMatrix, [], Cost, ecl, ConLev).

ham_path_g_c(Start, End, Succ, CostMatrix, Cost, ConLev) :-
        ham_path_offset_body(Start, End, Succ, 0, CostMatrix, [], Cost, gc, ConLev).


ham_path_offset(Start, End, Succ, Off, CostMatrix, Cost) :-
        ham_path_offset_body(Start, End, Succ, Off, CostMatrix, [], Cost, ecl, default).

ham_path_offset_g(Start, End, Succ, Off, CostMatrix, Cost) :-
        ham_path_offset_body(Start, End, Succ, Off, CostMatrix, [], Cost, gc, default).

ham_path_offset_c(Start, End, Succ, Off, CostMatrix, Cost, ConLev) :-
        ham_path_offset_body(Start, End, Succ, Off, CostMatrix, [], Cost, ecl, ConLev).

ham_path_offset_g_c(Start, End, Succ, Off, CostMatrix, Cost, ConLev) :-
        ham_path_offset_body(Start, End, Succ, Off, CostMatrix, [], Cost, gc, ConLev).


ham_path(Start, End, Succ, CostMatrix, ArcCosts, Cost) :-
        ham_path_offset_body(Start, End, Succ, 0, CostMatrix, ArcCosts, Cost, ecl, default).

ham_path_g(Start, End, Succ, CostMatrix, ArcCosts, Cost) :-
        ham_path_offset_body(Start, End, Succ, 0, CostMatrix, ArcCosts, Cost, gc, default).

ham_path_c(Start, End, Succ, CostMatrix, ArcCosts, Cost, ConLev) :-
        ham_path_offset_body(Start, End, Succ, 0, CostMatrix, ArcCosts, Cost, ecl, ConLev).

ham_path_g_c(Start, End, Succ, CostMatrix, ArcCosts, Cost, ConLev) :-
        ham_path_offset_body(Start, End, Succ, 0, CostMatrix, ArcCosts, Cost, gc, ConLev).


ham_path_offset(Start, End, Succ, Off, CostMatrix, ArcCosts, Cost) :-
        ham_path_offset_body(Start, End, Succ, Off, CostMatrix, ArcCosts, Cost, ecl, default).

ham_path_offset_g(Start, End, Succ, Off, CostMatrix, ArcCosts, Cost) :-
        ham_path_offset_body(Start, End, Succ, Off, CostMatrix, ArcCosts, Cost, gc, default).

ham_path_offset_c(Start, End, Succ, Off, CostMatrix, ArcCosts, Cost, ConLev) :-
        ham_path_offset_body(Start, End, Succ, Off, CostMatrix, ArcCosts, Cost, ecl, ConLev).

ham_path_offset_g_c(Start, End, Succ, Off, CostMatrix, ArcCosts, Cost, ConLev) :-
        ham_path_offset_body(Start, End, Succ, Off, CostMatrix, ArcCosts, Cost, gc, ConLev).


ham_path_offset_body(Start, End, Succ, Off0, CostMatrix, ArcCosts, Cost, IdxType, ConLev) :-
        check_integer(Off0),
        check_collection_to_list(Succ, SList),
        length(SList, N), 
        (IdxType == ecl ->
            Off is Off0 + 1, Min = Off, Max is N + Off 
        ;
            Off is Off0, Min = Off, Max is N + Off 
        ),
        get_prob_handle_nvars(H, NV0),
        process_domain_vars([Start,End|SList], [Min..Max], H, NV0,NV1, 
                            [],OldVs, [GStart,GEnd|GSs]),
        SArr =.. [[]|GSs],
        ham_path_offset_events1(GStart, GEnd, SArr, Off, CostMatrix,
                                ArcCosts, Cost, ConLev, H, NV1,NV, Ev,EvT),
        !,
        assign_domain([Min..Max], H, NV1, OldVs),
        % NV1 because assign_domain updated nvars to NV1
        update_gecode_with_default_newvars(H, NV1, NV),
        post_new_event_with_aux(Ev, EvT, H).
ham_path_offset_body(Start, End, Succ, Off0, CostMatrix, ArcCosts,
                     Cost, IdxType, _ConLev) :-
        get_bip_error(E),
        (IdxType = ecl ->
            (ArcCosts == [] ->
                error(E, ham_path_offset(Start, End, Succ, Off0, ArcCosts, Cost))
            ;
                error(E, ham_path_offset(Start, End, Succ, Off0, CostMatrix, ArcCosts, Cost))
            )
        ;
            (ArcCosts == [] ->
                error(E, ham_path_offset_g(Start, End, Succ, Off0, ArcCosts, Cost))
            ;
                error(E, ham_path_offset_g(Start, End, Succ, Off0, CostMatrix, ArcCosts, Cost))
            )
        ).

ham_path_offset_events1(GStart, GEnd, SArr, Off, CostMatrix, ArcCosts,
                        Cost, ConLev, H, NV1,NV, Ev,EvT) :-
        check_collection_to_list(flatten(CostMatrix), CMList),
        CMArr =.. [[]|CMList],
        ec_to_gecode_var1(Cost, H, NV1,NV2, GCost),
        ( ArcCosts == [] ->
            NV = NV2,
            GAC = []
        ;
            check_collection_to_list(ArcCosts, ACList),
            ec_to_gecode_varlist1(ACList, H, NV2,NV, GACList, _),
            GAC =.. [[]|GACList]
        ),
        Ev = [post_ham_path_cost(ConLev,GStart,GEnd,SArr,CMArr,GAC,GCost,Off)|EvT].

% The common code for posting a aux. ham_path (with cost) constraint in a
% constraint expression, we cannot treat the Succ variables differently here
ham_path_aux_events1(Start, End, Succ, Off, CostMat, ArcCosts, Cost, ConLev, H, N0,N, Auxs0,AuxsT) :-
        check_collection_to_list(Succ, SuccL),
        ec_to_gecode_varlist1([Start,End|SuccL], H, N0,N1, [GStart,GEnd|GSuccL], _),
        SArr =.. [[]|GSuccL],
        ham_path_offset_events1(GStart, GEnd, SArr, Off, CostMat, ArcCosts, Cost, ConLev, H, N1,N, Auxs0,AuxsT). 


sequence(Lo, Hi, K, Vars, Values) :-
        sequence_c(Lo, Hi, K, Vars, Values, default).

sequence_c(Lo, Hi, K, Vars, Values, ConLev) :-
        check_integer(Lo), check_integer(Hi), 
        Hi >= Lo, Lo >= 0,
        check_integer(K), K > 0,
        check_collection_to_list(Vars, VarList),
        check_collection_to_list(Values, ValList),
        (foreach(V, ValList) do check_integer(V)),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1(VarList, H, NV0,NV, GVars, _), !,
        update_gecode_with_default_newvars(H, NV0, NV),
        VarArr =.. [[]|GVars],
        ValArr =.. [[]|ValList],
        K =< arity(VarArr),
        post_new_event(post_sequence(ConLev, Lo, Hi, K, VarArr, ValArr), H).
sequence_c(Lo, Hi, K, Vars, Values, _ConLev) :-
        get_bip_error(E),
        error(E,sequence(Lo, Hi, K, Vars, Values)).


sequence(Lo, Hi, K, ZeroOnes) :-
        sequence_c(Lo, Hi, K, ZeroOnes, default).

sequence_c(Lo, Hi, K, ZeroOneVars, ConLev) :-
        check_integer(Lo), check_integer(Hi), 
        Hi >= Lo, Lo >= 0,
        check_integer(K), K > 0,
        check_collection_to_list(ZeroOneVars, VarList),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_boolvarlist1(VarList, H, NewVs, NBVs, NV0,NV, GVList, HasVar), !, 
        update_gecode_with_boolvars(NewVs, NBVs, NV0, NV, HasVar, H),
        VarArr =.. [[]|GVList],
        K =< arity(VarArr),
        post_new_event(post_sequence_01(ConLev, Lo, Hi, K, VarArr), H).
sequence_c(Lo, Hi, K, ZeroOneVars, _ConLev) :-
        get_bip_error(E),
        error(E,sequence_c(Lo, Hi, K, ZeroOneVars)).


cumulatives_min(Starts, Durations, Usages, UsedMachines, MachineMin) :-
        cumulatives_c(Starts, Durations, Usages, UsedMachines, MachineMin, 0,
                      ecl, default).

cumulatives_min_g(Starts, Durations, Usages, UsedMachines, MachineMin) :-
        cumulatives_c(Starts, Durations, Usages, UsedMachines, MachineMin, 0,
                      gc, default).

cumulatives(Starts, Durations, Usages, UsedMachines, MachineLimits) :-
        cumulatives_c(Starts, Durations, Usages, UsedMachines, MachineLimits,
                      1, ecl, default).

cumulatives_g(Starts, Durations, Usages, UsedMachines, MachineLimits) :-
        cumulatives_c(Starts, Durations, Usages, UsedMachines, MachineLimits,
                      1, gc, default).

cumulatives_min_c(Starts, Durations, Usages, UsedMachines, MachineMin, ConLev) :-
        cumulatives_c(Starts, Durations, Usages, UsedMachines, MachineMin, 0, ecl, ConLev).

cumulatives_min_g_c(Starts, Durations, Usages, UsedMachines, MachineMin, ConLev) :-
        cumulatives_c(Starts, Durations, Usages, UsedMachines, MachineMin, 0, gc, ConLev).

cumulatives_c(Starts, Durations, Usages, UsedMachines, MachineMin, ConLev) :-
        cumulatives_c(Starts, Durations, Usages, UsedMachines, MachineMin, 1, ecl, ConLev).

cumulatives_g_c(Starts, Durations, Usages, UsedMachines, MachineMin, ConLev) :-
        cumulatives_c(Starts, Durations, Usages, UsedMachines, MachineMin, 1, gc, ConLev).

cumulatives_c(Starts, Durations, Usages, UsedMachines, Limits, AtMost, IndexType, ConLev) :-
        check_collection_to_list(Starts, StartsL),
        check_collection_to_list(Durations, DurationsL),
        check_collection_to_list(Usages, UsagesL),
        check_collection_to_list(UsedMachines, UsedL),
        get_prob_handle_nvars(H, NV0),
        H = gfd_prob{nvars:NV0},
        ec_to_gecode_multivarlists1([StartsL,DurationsL,UsagesL,UsedL], 
                                    H, NV0,NV1,
                                   [GStartsL,GDurationsL,GUsagesL,GUsedL]),
        StartsArr =.. [[]|GStartsL],
        arity(StartsArr, N),
        DurationsArr =.. [[]|GDurationsL],
        arity(DurationsArr,N),
        UsagesArr =.. [[]|GUsagesL],
        arity(UsagesArr,N),
        UsedArr =.. [[]|GUsedL],
        arity(UsedArr,N),
        functor(EndsArr, [], N),
        create_task_end_times1(StartsArr, DurationsArr, H, NV1,NV, Evs, EvsT0, EndsArr),
        check_collection_to_list(Limits, LimitsL),
	(IndexType == ecl ->
            LimitsArr =.. [[],0|LimitsL], % dummy limit for 0'th machine 
	    First = 1
        ;
            LimitsArr =.. [[]|LimitsL],
	    First = 0
	), !,

        % must have new vars as Ends' vars are created by constraint
        gfd_default_interval(Min, Max),
        do_update_newvars_with_domain_interval(H, NV, Min, Max),
        
        Last is arity(LimitsArr) - 1,
        EvsT0 = [post_cumulatives(ConLev, StartsArr, DurationsArr, EndsArr,
                                  UsagesArr, UsedArr, LimitsArr, AtMost), 
                 post_interval(UsedArr,First,Last) |
                 EvsT],
        post_new_event_with_aux(Evs, EvsT, H).
cumulatives_c(Starts, Durations, Usages, UsedMachines, Limits, AtMost, IndexType, _ConLev) :-
        get_bip_error(E),
        (IndexType = ecl ->
            error(E,cumulatives(Starts, Durations, Usages,
                                UsedMachines, Limits, AtMost))
        ;
            error(E,cumulatives_g(Starts, Durations, Usages,
                                  UsedMachines, Limits, AtMost))
        ).


cumulative(Starts, Durations, Usages, Limit) :-
        get_prob_handle(H),
        cumulative_body(Starts, Durations, Usages, Limit, H, []).

cumulative_optional(Starts, Durations, Usages, Limit, Scheduled) :-
        check_collection_to_list(Scheduled, BList),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_boolvarlist1(BList, H, NewVs, NBVs, NV0,NV, GBs, HasVar),
        !,
        update_gecode_with_boolvars(NewVs, NBVs, NV0, NV, HasVar, H),
        cumulative_body(Starts, Durations, Usages, Limit, H, GBs).
cumulative_optional(Starts, Durations, Usages, Limit, Scheduled) :-
        get_bip_error(E),
        error(E, cumulative_optional(Starts, Durations, Usages, Limit, Scheduled)).
        
cumulative_body(Starts, Durations, Usages, Limit, H, GBs) :-
        check_collection_to_list(Starts, StartsL),
        check_collection_to_list(Durations, DurationsL),
        H = gfd_prob{nvars:NV0},
        ec_to_gecode_varlist1([Limit|StartsL], H, NV0,NV1, [GLimit|GStartsL], _),
        ec_to_gecode_varlist1(DurationsL, H, NV1,NV2, GDsL, DHasVar),
        SsArr =.. [[]|GStartsL],
        arity(SsArr, N),
        DsArr =.. [[]|GDsL],
        arity(DsArr,N),
        % usages must be integers
        check_collection_to_list(Usages, UsagesL),
        UsArr =.. [[]|UsagesL],
        arity(UsArr,N),

        ( GBs == [] ->
            BsArr = []
        ;
            BsArr =.. [[]|GBs],
            arity(BsArr, N)
        ), !,
        ( var(DHasVar) ->
            update_gecode_with_default_newvars(H, NV0, NV2),
            post_new_event(post_cumulative(SsArr, DsArr, UsArr, GLimit, BsArr), H) 
        ;
            dim(EsArr, [N]),
            create_task_end_times1(SsArr, DsArr, H, NV2,NV, Auxs0,AuxsT, EsArr),
            update_gecode_with_default_newvars(H, NV0, NV),
            post_new_event_with_aux([post_cumulativeflex(SsArr, DsArr, EsArr, UsArr, GLimit, BsArr)|Auxs0], AuxsT, H)                                     
        ).
cumulative_body(Starts, Durations, Usages, Limit, H, GBs) :-
        get_bip_error(E),
        % too difficult to reconstruct original cumulative call, 
        % just raise error with cumulative_body
        error(E, cumulative_body(Starts, Durations, Usages, Limit, H, GBs)).


disjunctive(Starts, Durations) :-
        get_prob_handle(H),
        disjunctive_body(Starts, Durations, H, []).

disjunctive_optional(Starts, Durations, Scheduled) :-
        get_prob_handle_nvars(H, NV0),
        check_collection_to_list(Scheduled, BList),
        ec_to_gecode_boolvarlist1(BList, H, NewVs, NBVs, NV0,NV, GBs, HasVar),
        !,
        update_gecode_with_boolvars(NewVs, NBVs, NV0, NV, HasVar, H),
        disjunctive_body(Starts, Durations, H, GBs).
disjunctive_optional(Starts, Durations, Scheduled) :-
        get_bip_error(E),
        error(E, disjunctive_optional(Starts, Durations, Scheduled)).

disjunctive_body(Starts, Durations, H, GBs) :-
        check_collection_to_list(Starts, SList),
        check_collection_to_list(Durations, DList),
        H = gfd_prob{nvars:NV0},
        ec_to_gecode_varlist1(SList, H, NV0,NV1, GSs, _),
        ec_to_gecode_varlist1(DList, H, NV1,NV2, GDs, DHasVar),
        SsArr =.. [[]|GSs],
        DsArr =.. [[]|GDs],
        arity(SsArr, N),
        arity(DsArr, N),  % same size
        ( GBs == [] ->
            BsArr = []
        ;
            BsArr =.. [[]|GBs],
            arity(BsArr, N)
        ), !,
        ( var(DHasVar) ->
            update_gecode_with_default_newvars(H, NV0, NV2),
            post_new_event(post_disj(SsArr,DsArr,BsArr), H)
        ;
            functor(EsArr, [], N),
            create_task_end_times1(SsArr, DsArr, H, NV2,NV, Auxs0,AuxsT, EsArr),
            update_gecode_with_default_newvars(H, NV0, NV),
            post_new_event_with_aux([post_disjflex(SsArr,DsArr,EsArr,BsArr)|Auxs0], AuxsT, H)
        ).
disjunctive_body(Starts, Durations, H, GBs) :-
        get_bip_error(E),
        error(E, disjunctive_body(Starts, Durations, H, GBs)).


create_task_end_times1(SsArr, DsArr, H, NV0,NV, Auxs0,AuxsT, EsArr) :-
        ( foreacharg(GS, SsArr), foreacharg(GD, DsArr),
          param(H),
          fromto(NV0, NV1,NV2, NV), fromto(Auxs0, Auxs1,Auxs2, AuxsT),
          foreacharg(GE, EsArr)
        do
%            ec_to_gecode_var1(E, H, NV1,NV2, GE),
            % EsArr are all new variables
            ( integer(GS), integer(GD) ->
                % don't post extra constraint if End could be determined
                GE is GS + GD,
                NV1 = NV2,
                Auxs1 = Auxs2
            ;
                new_gfdvar(_E, H, NV1,NV2, GE),
                Auxs1 = [post_sum(default, [](GS,GD), (#=), GE)|Auxs2]
            )
        ).


bool_channeling(V, Bools, Min) :-
        bool_channeling_c(V, Bools, Min, default).

bool_channeling_c(V, Bools, Min, ConLev) :-
        check_integer(Min),
        check_collection_to_list(Bools, BList),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_boolvarlist1(BList, H, NewVs, NBVs, NV0,NV1, GBList, HasVar),
        !,
        update_gecode_with_boolvars(NewVs, NBVs, NV0, NV1, HasVar, H),
        GBArr =.. [[]|GBList],
        Max is Min+arity(GBArr),
        ( var(V) ->
            ec_to_gecode_oldvar1(V, H, NV1,NV, [], OldGV, GV),
            assign_domain_interval1(H, NV1, NV, OldGV, Min, Max)
        ; integer(V) ->
            V >= Min, 
            V =< Max,
            GV = V
        ;
            error(5, bool_channeling(V, Bools, Min))
        ),
        post_new_event(post_boolchannel(ConLev,GV,GBArr,Min), H). 
bool_channeling_c(V, Bools, Min, _ConLev) :-
        get_bip_error(E),
        error(E, bool_channeling(V, Bools, Min)).


:- export struct(rect(x,y,w,h,b)).

disjoint2(Recs) :-
        check_collection_to_list(Recs, RecLs),
        ( foreach(rect{x:X,y:Y,w:Width,h:Height}, RecLs),
          foreach(X, Xs), foreach(Y, Ys), 
          foreach(Width, Ws), foreach(Height, Hs)
        do 
            true
        ),
        get_prob_handle(H),
        disjoint2_body(Xs, Ws, Ys, Hs, [], H).


disjoint2_optional(Recs) :-
        check_collection_to_list(Recs, RecLs),
        ( foreach(rect{x:X,y:Y,w:Width,h:Height,b:O}, RecLs),
          foreach(X, Xs), foreach(Y, Ys), foreach(O, Os), 
          foreach(Width, Ws), foreach(Height, Hs)
        do 
            true
        ),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_boolvarlist1(Os, H, NewVs, NBVs, NV0, NV, GOs, HasVar),
        !,
        update_gecode_with_boolvars(NewVs, NBVs, NV0, NV, HasVar, H),
        disjoint2_body(Xs, Ws, Ys, Hs, GOs, H).
disjoint2_optional(Recs) :-
        get_bip_error(E),
        error(E, disjoint2_optional(Recs)).

disjoint2_body(Xs, Ws, Ys, Hs, GOs, H) :-
        H = gfd_prob{nvars:N0}, 
        ec_to_gecode_varlist1(Xs, H, N0,N1, GXs, _),
        ec_to_gecode_varlist1(Ys, H, N1,N2, GYs, _),
        ec_to_gecode_varlist1(Hs, H, N2,N3, GHs, HasVar),
        ec_to_gecode_varlist1(Ws, H, N3,N4, GWs, HasVar),
        WArr =.. [[]|GWs],
        HArr =.. [[]|GHs],
        XArr =.. [[]|GXs],
        YArr =.. [[]|GYs],
        OArr =.. [[]|GOs],
        !,
        ( HasVar == 1 ->
            arity(WArr, NR),
            dim(XEArr, [NR]),
            dim(YEArr, [NR]),
            create_task_end_times1(XArr, WArr, H, N4, N5, Auxs0,Auxs1, XEArr),
            create_task_end_times1(YArr, HArr, H, N5, N6, Auxs1,AuxsT, YEArr),
            update_gecode_with_default_newvars(H, N0, N6),
            post_new_event_with_aux([post_disjointflex2(XArr,WArr,YArr,HArr,OArr,XEArr,YEArr)|
                                     Auxs0], AuxsT, H)
        ;
            update_gecode_with_default_newvars(H, N0, N4),
            post_new_event(post_disjoint2(XArr,WArr,YArr,HArr,OArr), H)
        ).
disjoint2_body(Xs, Ws, Ys, Hs, GOs, H) :-
        get_bip_error(E),
        error(E, disjoint2_body(Xs, Ws, Ys, Hs, GOs, H)).


inverse(XL, YL) :-
        inverse_body(XL, YL, ecl, default).

inverse_c(XL, YL, ConLev) :-
        inverse_body(XL, YL, ecl, ConLev).

inverse_g(XL, YL) :-
        inverse_body(XL, YL, gc, default).

inverse_g_c(XL, YL, ConLev) :-
        inverse_body(XL, YL, gc, ConLev).

inverse_body(Vs1, Vs2, IndexType, ConLev) :-
        ( var(Vs1) -> 
            Vs1 = Vars1
        ;
            check_collection_to_list(Vs1, Vars1)
        ),
        ( var(Vs2) -> 
            nonvar(Vs1),
            Vs2 = Vars2
        ;
            check_collection_to_list(Vs2, Vars2)
        ),
	( foreach(_,Vars1), foreach(_,Vars2) do true ),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1(Vars1, H, NV0,NV1, GVars1, _),
        ec_to_gecode_varlist1(Vars2, H, NV1,NV, GVars2, _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        % add 0 to match ECLiPSe index (starting from 1)
        (IndexType == ecl ->
            Arr1 =.. [[],0|GVars1],
            Arr2 =.. [[],0|GVars2]
        ;
            Arr1 =.. [[]|GVars1],
            Arr2 =.. [[]|GVars2]
        ),
        post_new_event(post_inverse(ConLev, Arr1, Arr2), H).
inverse_body(Vs1, Vs2, IndexType, _ConLev) :-
        get_bip_error(E),
        (IndexType = ecl ->
            error(E, inverse(Vs1, Vs2))
        ;
            error(E, inverse_g(Vs1, Vs2))
        ).


inverse(XL, XOff, YL, YOff) :-
        inverse_body(XL, XOff, YL, YOff, ecl, default).

inverse_c(XL, XOff, YL, YOff, ConLev) :-
        inverse_body(XL, XOff, YL, YOff, ecl, ConLev).

inverse_g(XL, XOff, YL, YOff) :-
        inverse_body(XL, XOff, YL, YOff, gc, default).

inverse_g_c(XL, XOff, YL, YOff, ConLev) :-
        inverse_body(XL, XOff, YL, YOff, gc, ConLev).

inverse_body(Vs1, Off1, Vs2, Off2, IndexType, ConLev) :-
        check_integer(Off1),
        check_integer(Off2),
	(IndexType == ecl ->
	    GOff1 is Off1 - 1,
            GOff2 is Off2 - 1
	;
	    GOff1 is Off1,
	    GOff2 is Off2
	),
        ( var(Vs1) -> 
            Vs1 = Vars1
        ;
            check_collection_to_list(Vs1, Vars1)
        ),
        ( var(Vs2) -> 
            nonvar(Vs1),
            Vs2 = Vars2
        ;
            check_collection_to_list(Vs2, Vars2)
        ),
	( foreach(_,Vars1), foreach(_,Vars2) do true ),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1(Vars1, H, NV0,NV1, GVars1, _),
        ec_to_gecode_varlist1(Vars2, H, NV1,NV, GVars2, _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        Arr1 =.. [[]|GVars1],
        Arr2 =.. [[]|GVars2],
        post_new_event(post_inverse_offset(ConLev, Arr1, GOff1, Arr2, GOff2), H).
inverse_body(Vs1, Off1, Vs2, Off2, IndexType, _ConLev) :-
        get_bip_error(E),
        (IndexType = ecl ->
            error(E, inverse(Vs1, Off1, Vs2, Off2))
        ;
            error(E, inverse_g(Vs1, Off1, Vs2, Off2))
        ).


precede(S, T, Vars) :-
        precede_c(S, T, Vars, default).

precede_c(S, T, Vars, _ConLev) :-  % ConLev ignored; constraint is gac
        check_integer(S),
        check_integer(T),
        check_collection_to_list(Vars, VList),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1(VList, H, NV0, NV, GVars, _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        GArr =.. [[]|GVars],
        post_new_event(post_precede(S,T,GArr), H).
precede_c(S, T, Vars, _ConLev) :-
        get_bip_error(E),
        error(E, precede(S, T, Vars)).


precede(Vals, Vars) :-
        check_collection_to_list(Vals, ValsL),
        ValsArr =.. [[]|ValsL],
        check_collection_to_list(Vars, VarsL),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1(VarsL, H, NV0,NV, GVarsL, _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        GVarsArr =.. [[]|GVarsL],
        post_new_event(post_precede_chain(ValsArr,GVarsArr), H).
precede(Vals, Vars) :-
        get_bip_error(E),
        error(E, precede(Vals, Vars)).
        

min(Xs, Min) :-
        minlist_c(Xs, Min, default).

minlist_c(Xs, Min, ConLev) :-
        check_collection_to_list(Xs, XLs),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([Min|XLs], H, NV0,NV, [GMin|GLs], _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        GArray =.. [[]|GLs],
        post_new_event(post_minlist(ConLev, GMin, GArray), H).
minlist_c(Xs, Min, _ConLev) :-
        get_bip_error(E),
        error(E, min(Xs, Min)).


max(Xs, Max) :-
        maxlist_c(Xs, Max, default).

maxlist_c(Xs, Max, ConLev) :-
        check_collection_to_list(Xs, XLs),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([Max|XLs], H, NV0,NV, [GMax|GLs], _),
        GArray =.. [[]|GLs],
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_maxlist(ConLev, GMax, GArray), H).
maxlist_c(Xs, Max, _ConLev) :-
        get_bip_error(E),
        error(E, max(Xs, Max)).


sumlist(Xs, Sum) :-
        sum_c(Xs, Sum, default).

sum(Xs, Sum) :-
        sum_c(Xs, Sum, default).

sum_c(Xs, Sum, ConLev) :-
        nonvar(Xs),
        ( Xs = Cs * Vs ->
            scalar_product_c(Cs, Vs, (#=), Sum, ConLev)
        ;
            sum_c(Xs, (#=), Sum, ConLev)
        ).


sum(Xs, Rel, Sum) :-
        sum_c(Xs, Rel, Sum, default).

sum_c(Xs, Rel, Sum, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        linsum_event1(Xs, Rel, Sum, ConLev, H, NV0,NV, Events,EventsT),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event_with_aux(Events, EventsT, H).
sum_c(Xs, Rel, Sum, _ConLev) :-
        get_bip_error(E),
        error(E, sum(Xs, Rel, Sum)).


linsum_event1(Xs, Rel0, Sum, ConLev, H, NV0,NV, Events,EventsT) :-
        linsum_body1(Xs, Rel0, Sum, GArray, Rel, GSum, H, NV0,NV),
        Events = [post_sum(ConLev, GArray, Rel, GSum)|EventsT].

linsum_body1(Xs, Rel0, Sum, GArray, Rel, GSum, H, NV0,NV) :-
	is_valid_rel_op(Rel0, Rel),
        check_collection_to_list(Xs, XLs),
        ec_to_gecode_varlist1([Sum|XLs], H, NV0,NV, [GSum|GLs], _),
        GArray =.. [[]|GLs].


sum(Xs, Rel, Sum, Bool) :-
        sum_reif_c(Xs, Rel, Sum, Bool, default).

sum_reif_c(Xs, Rel, Sum, Bool, ConLev) :-
	get_prob_handle_nvars(H, NV0),
	linsum_reif_event1(Xs, Rel, Sum, Bool, ConLev, H, NV0,NV, [],
                           Bs, Event, _GBool),
        !,
	gfd_default_interval(Min, Max),
	update_vars_for_gecode(NV0, NV, Bs, H, Min, Max),
        post_new_event_with_aux([Event|EventsT], EventsT, H).
sum_reif_c(Xs, Rel, Sum, Bool, _ConLev) :-
        get_bip_error(E),
        error(E, sum(Xs, Rel, Sum, Bool)).

linsum_reif_event1(Xs, Rel0, Sum, Bool, ConLev, H, NV0,NV, Bs0,Bs, Event, GBool) :-
        linsum_body1(Xs, Rel0, Sum, GArray, Rel, GSum, H, NV0,NV1),
        ( var(Bool) ->
            ec_to_gecode_var1(Bool, H, NV1,NV, GBool),
            Bs = [Bool|Bs0]
        ;
            GBool = Bool,
            NV1 = NV,
            Bs = Bs0
        ;
            fail
        ),
        Event = post_sum_reif(ConLev, GArray, Rel, GSum, GBool).


scalar_product(Cs, Xs, Rel, P) :-
        scalar_product_c(Cs, Xs, Rel, P, default).

scalar_product_c(Cs, Xs, Rel, P, ConLev) :-
	get_prob_handle_nvars(H, NV0),
	scalar_product_event1(Cs, Xs, Rel, P, ConLev, H, NV0,NV, Events,EventsT),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event_with_aux(Events, EventsT, H).
scalar_product_c(Cs, Xs, Rel, P, _ConLev) :-
        get_bip_error(E),
        error(E, scalar_product_c(Cs, Xs, Rel, P)).

scalar_product_event1(Cs, Xs, Rel0, P, ConLev, H, NV0,NV, Events,EventsT) :-
	scalar_product_body1(Cs, Xs, Rel0, P, CArray, GArray, GP, Rel, H, NV0,NV),
	Events = [post_lin(ConLev, GArray, CArray, Rel, GP)|EventsT].
		  
scalar_product_body1(Cs, Xs, Rel0, P, CArray, GArray, GP, Rel, H, NV0,NV) :-
	is_valid_rel_op(Rel0, Rel),
        check_collection_to_list(Xs, XLs),
        check_collection_to_list(Cs, CLs),
        ec_to_gecode_varlist1([P|XLs], H, NV0,NV, [GP|GLs], _),
        GArray =.. [[]|GLs],
	CArray =.. [[]|CLs],
	arity(GArray) =:= arity(CArray),
	(foreach(C, CLs) do check_integer(C)).


scalar_product(Cs, Xs, Rel, P, Bool) :-
        scalar_product_reif_c(Cs, Xs, Rel, P, Bool, default).


scalar_product_reif_c(Cs, Xs, Rel, P, Bool, ConLev) :-
	get_prob_handle_nvars(H, NV0),
	scalar_product_reif_event1(Cs, Xs, Rel, P, Bool, ConLev, H,
                                   NV0,NV, [],Bs, Event, _GBool),
        !,
	gfd_default_interval(Min, Max),
	update_vars_for_gecode(NV0, NV, Bs, H, Min, Max),
        post_new_event_with_aux([Event|EventsT], EventsT, H).
scalar_product_reif_c(Cs, Xs, Rel, P, Bool, _ConLev) :-
        get_bip_error(E),
        error(E, scalar_product(Cs, Xs, Rel, P, Bool)).


scalar_product_reif_event1(Cs, Xs, Rel0, P, Bool, ConLev, H, NV0,NV, Bs0,Bs, Event, GBool) :-
	scalar_product_body1(Cs, Xs, Rel0, P, CArray, GArray, GP, Rel,
                             H, NV0,NV1),
        ( var(Bool) ->
            ec_to_gecode_var1(Bool, H, NV1,NV, GBool),
            Bs = [Bool|Bs0]
        ; integer(Bool) ->
            Bool >= 0,
            Bool =< 1,
            GBool = Bool,
            NV1 = NV,
            Bs = Bs0
        ;
            set_bip_error(5)
        ),
	Event = post_lin_reif(ConLev, GArray, CArray, Rel, GP, GBool).


% X #= Y + C => X - Y #= C
ac_eq(X, Y, C) :-
        integer(C),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([X,Y], H, NV0,NV, GXY, _),
        GArray =.. [[]|GXY],
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_lin(gfd_gac, GArray, [](1,-1), (#=), C), H).



isqrt_c(X, Y, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([X,Y], H, NV0,NV, [GX,GY], _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_sqrt(ConLev, GY,GX), H).
isqrt_c(X, Y, _ConLev) :-
        get_bip_error(E),
        error(E, isqrt(X, Y)).

sqr(X, Y) :-
        sqr_c(X, Y, default).

sqr_c(X, Y, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([X,Y], H, NV0,NV, [GX,GY], _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_sq(ConLev, GY,GX), H).
sqr_c(X, Y, _ConLev) :-
        get_bip_error(E),
        error(E, sqr(X, Y)).


abs_c(X, Y, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([X,Y], H, NV0,NV, [GX,GY], _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_abs(ConLev, GY,GX), H).
abs_c(X, Y, _ConLev) :-
        get_bip_error(E),
        error(E, abs(X, Y)).


div_c(X, Y, Z, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([X,Y,Z], H, NV0,NV, [GX,GY,GZ], _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_div(ConLev, GZ, GX,GY), H).
div_c(X, Y, Z, _ConLev) :-
        get_bip_error(E),
        error(E, div(X, Y, Z)).


divmod(X, Y, Q, M) :-
	divmod_c(X, Y, Q, M, default).

divmod_c(X, Y, Q, M, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([X,Y,Q,M], H, NV0,NV, [GX,GY,GQ,GM], _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_divmod(ConLev, GX, GY, GQ, GM), H).
divmod_c(X, Y, Q, M, _ConLev) :-
        get_bip_error(E),
        error(E, divmod_c(X, Y, Q, M)).


mult_c(X, Y, Z, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([X,Y,Z], H, NV0,NV, [GX,GY,GZ], _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_mult(ConLev, GZ, GX,GY), H).
mult_c(X, Y, Z, _ConLev) :-
        get_bip_error(E),
        error(E, mult(X, Y, Z)).


mod_c(X, Y, Z, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([X,Y,Z], H, NV0,NV, [GX,GY,GZ], _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_mod(ConLev, GZ, GX,GY), H).
mod_c(X, Y, Z, _ConLev) :-
        get_bip_error(E),
        error(E, rem(X, Y, Z)).



min_c(X, Y, Z, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([X,Y,Z], H, NV0,NV, [GX,GY,GZ], _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_min2(ConLev, GZ, GX,GY), H).
min_c(X, Y, Z, _ConLev) :-
        get_bip_error(E),
        error(E, min(X, Y, Z)).


max_c(X, Y, Z, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([X,Y,Z], H, NV0,NV, [GX,GY,GZ], _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_max2(ConLev, GZ, GX,GY), H).
max_c(X, Y, Z, _ConLev) :-
        get_bip_error(E),
        error(E, max(X, Y, Z)).


ordered(Order, Xs) :-
        ordered_c(Order, Xs, default).

ordered_c(Order0, Xs0, ConLev) :-
        is_valid_rel_op(Order0, Order),
        get_prob_handle_nvars(H, NV0),
        check_collection_to_list(flatten(Xs0), Xs),
        ec_to_gecode_varlist1(Xs, H, NV0,NV, GXs, _),
        XArr =.. [[]|GXs],
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_ordered(ConLev, XArr, Order), H).
ordered_c(Order0, Xs0, _ConLev) :-
        get_bip_error(E),
        error(E, ordered(Order0, Xs0)).


all_le(X, Y) :-
        collection_rel(X, (#=<), Y, default).

all_le_c(X, Y, ConLev) :-
        collection_rel(X, (#=<), Y, ConLev).

all_lt(X, Y) :-
        collection_rel(X, (#<), Y, default).

all_lt_c(X, Y, ConLev) :-
        collection_rel(X, (#<), Y, ConLev).

all_ge(X, Y) :-
        collection_rel(X, (#>=), Y, default).

all_ge_c(X, Y, ConLev) :-
        collection_rel(X, (#>=), Y, ConLev).

all_gt(X, Y) :-
        collection_rel(X, (#>), Y, default).

all_gt_c(X, Y, ConLev) :-
        collection_rel(X, (#>), Y, ConLev).

all_ne(X, Y) :-
        collection_rel(X, (#\=), Y, default).

all_ne_c(X, Y, ConLev) :-
        collection_rel(X, (#\=), Y, ConLev).

all_eq(X, Y) :-
        collection_rel(X, (#=), Y, default).

all_eq_c(X, Y, ConLev) :-
        collection_rel(X, (#=), Y, ConLev).


rel_c(X, Rel, Y, ConLev) :-
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1([X,Y], H, NV0,NV, [GX,GY], _),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_rel(ConLev, GX, Rel, GY), H).
rel_c(X, Rel, Y, _ConLev) :-
        get_bip_error(E),
        error(E, rel_c(X, Rel, Y)).


collection_rel(Xs, Rel, Y, ConLev) :-
        check_collection_to_list(Xs, XList),
        (XList == [] ->
            !,
            true
        ;
            get_prob_handle_nvars(H, NV0),
            ec_to_gecode_varlist1([Y|XList], H, NV0,NV, [GY|GXList], _),
            GXArr =.. [[]|GXList],
            !,
            update_gecode_with_default_newvars(H, NV0, NV),
            post_new_event(post_collection_rel(ConLev, GXArr, Rel, GY),H)
        ).
collection_rel(Xs, Rel, Y, _ConLev) :-
        get_bip_error(E),
        error(E, rel(Xs, Rel, Y)).

        
lex_le(Xs, Ys) :-
        lex_le_c(Xs, Ys, default).

lex_le_c(Xs, Ys, ConLev) :-
        lex_c(Xs, (#=<), Ys, ConLev).


lex_lt(Xs, Ys) :-
        lex_lt_c(Xs, Ys, default).

lex_lt_c(Xs, Ys, ConLev) :-
        lex_c(Xs, (#<), Ys, ConLev).


lex_ge(Xs, Ys) :-
        lex_ge_c(Xs, Ys, default).

lex_ge_c(Xs, Ys, ConLev) :-
        lex_c(Xs, (#>=), Ys, ConLev).


lex_gt(Xs, Ys) :-
        lex_gt_c(Xs, Ys, default).

lex_gt_c(Xs, Ys, ConLev) :-
        lex_c(Xs, (#>), Ys, ConLev).


lex_eq(Xs, Ys) :-
        lex_eq_c(Xs, Ys, default).

lex_eq_c(Xs, Ys, ConLev) :-
        lex_c(Xs, (#=), Ys, ConLev).


lex_ne(Xs, Ys) :-
        lex_ne_c(Xs, Ys, default).

lex_ne_c(Xs, Ys, ConLev) :-
        lex_c(Xs, (#\=), Ys, ConLev).


lex_c(Xs0, Rel0, Ys0, ConLev) :-
        is_valid_rel_op(Rel0, Rel),
        check_collection_to_list(flatten(Xs0), Xs),
        check_collection_to_list(flatten(Ys0), Ys),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1(Xs, H, NV0,NV1, GXs, _),
        XArr =.. [[]|GXs],
        ec_to_gecode_varlist1(Ys, H, NV1,NV, GYs, _),
        YArr =.. [[]|GYs],
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_new_event(post_lex_order(ConLev, XArr, Rel, YArr), H).
lex_c(Xs0, Rel0, Ys0, _ConLev) :-
        get_bip_error(E),
        error(E, lex(Xs0, Rel0, Ys0)).


bin_packing(Items, ItemSizes, BinLoads) :-
        bin_packing_body(Items, ItemSizes, BinLoads, ecl).

bin_packing_g(Items, ItemSizes, BinLoads) :-
        bin_packing_body(Items, ItemSizes, BinLoads, gc).

bin_packing_body(Items0, ItemSizes0, BinLoads0, IndexType) :-
        check_collection_to_list(flatten(Items0), Items),
        check_collection_to_list(flatten(ItemSizes0), ItemSizes),
        check_collection_to_list(flatten(BinLoads0), BinLoads),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_multivarlists1([Items,BinLoads], H, NV0,NV, [GIs,GLs]),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        post_bin_packing(GIs, ItemSizes, GLs, H, IndexType).
bin_packing_body(Items0, ItemSizes0, BinLoads0, IndexType) :-
        get_bip_error(E),
        (IndexType = ecl ->
            error(E, bin_packing(Items0, ItemSizes0, BinLoads0))
        ;
            error(E, bin_packing_g(Items0, ItemSizes0, BinLoads0))
        ).


post_bin_packing(GIs, ItemSizes, GLs, H, IndexType) :-
        (IndexType == ecl ->
            IArray =.. [[],0|GIs],
            SArray =.. [[],1|ItemSizes],
            arity(IArray) =:= arity(SArray),
            LArray =.. [[],1|GLs]
        ;
            IArray =.. [[]|GIs],
            SArray =.. [[]|ItemSizes],
            arity(IArray) =:= arity(SArray),
            LArray =.. [[]|GLs]
        ),
        post_new_event(post_bin_packing(IArray, SArray, LArray), H).


bin_packing(Items0,ItemSizes0,N,BinSize):-
        get_prob_handle_nvars(H, NV0),
        ( count(_, 1, N), 
          param(H), 
          fromto(NV0, NV1,NV2, NV3), foreach(GB, GBinLoads)
        do
            new_gfdvar(_B, H, NV1,NV2, GB)
        ),
        check_collection_to_list(flatten(Items0), Items),
        check_collection_to_list(flatten(ItemSizes0), ItemSizes),
        ec_to_gecode_varlist1(Items, H, NV3,NV, GItems, _),
        !,
        do_update_newvars_with_domain_interval(H, NV3, 0, BinSize),
        update_gecode_with_default_newvars(H, NV3, NV),
        post_bin_packing(GItems, ItemSizes, GBinLoads, H, ecl).
bin_packing(Items0,ItemSizes0,N,BinSize):-
        get_bip_error(E),
        error(E, bin_packing(Items0,ItemSizes0,N,BinSize)).


table(Vars, Table) :-
        table_c(Vars, Table, default, default).

table(Vars, Table, Emph) :-
        table_c(Vars, Table, Emph, default).

table_c(Vars, Table, ConLev) :-
        table_c(Vars, Table, default, ConLev).

table_c(Vars, Table, Emph, ConLev) :-
        check_table_emph(Emph),
        check_compound_to_list(Vars, VList),
        check_compound_to_list(Table, TList0),
        ( foreach(T0, TList0), param(TSize),
          foreach(T,TList)
         do
            check_nonvar(T0),
            ( T0 = [_|_] -> 
                T =.. [[]|T0]
            ; atomic(T0) -> 
                set_bip_error(5)
            ;
                T0 = T
            ),
            (arity(T, TSize) -> true ; set_bip_error(6))
        ),
        get_prob_handle_nvars(H, NV0),
        ( VList = [E|_], compound(E) ->
            /* multiple tuples of variables */
            ( foreach(Vs, VList), fromto(NV0, NV1,NV2, NV),
              param(H, TSize),
              foreach(GVArr, GVs)
            do
                check_collection_to_list(Vs, VsList),
                ec_to_gecode_varlist1(VsList, H, NV1,NV2, GVsList, _),
                GVArr =.. [[]|GVsList],
                ( arity(GVArr, TSize) -> true ; set_bip_error(6))
            )
        ;
            /* single tuple of variables */
            ec_to_gecode_varlist1(VList, H, NV0,NV, GVList, _),
            VArr =.. [[]|GVList],
            (arity(VArr, TSize) -> true ; set_bip_error(6)),
            GVs = [VArr]
        ), !,
        update_gecode_with_default_newvars(H, NV0, NV),
        g_create_tupleset_handle(TList, TSize, TsH),
        post_new_event(post_table(ConLev, GVs, TsH, TSize, Emph), H).
table_c(Vars, Table, Emph, _ConLev) :-
        get_bip_error(E),
        error(E, table(Vars, Emph, Table)).

check_table_emph(mem) ?- !.
check_table_emph(speed) ?- !.
check_table_emph(default) ?- !.
check_table_emph(_) :- set_bip_error(6).


% the order of the fields follows that used by Gecode's DFA::Transition
:- export struct(trans(f,l,t)).

extensional(Vars, Transitions, Start, Finals) :-
        extensional_c(Vars, Transitions, Start, Finals, default).

extensional_c(Vars, Transitions, Start, Finals, ConLev) :-
        check_integer(Start),
        check_collection_to_list(Finals, FList0),
        check_compound_to_list(Vars, VList),
        check_collection_to_list(Transitions, TList0),
        ( foreach(T, TList0), fromto(TList, TL1,TL2, TTail) do
            TL1 = [T|TL2],
            (T = trans{f:Fr,t:To,l:Sy},
              integer(Fr), integer(To), integer(Sy) ->
                check_nonnegative(Fr), check_nonnegative(To)
            ;
                set_bip_error(5)
            )

        ),
        TTail = [trans{f: -1, t:0,l:0}], % dummy transition to mark end
        TArr =.. [[]|TList],
        ( foreach(F, FList0), 
          fromto(FList, [F|T],T, [-1]) % as FList0 with extra -1 at end
        do 
            check_integer(F), check_nonnegative(F)
        ),
        FArr =.. [[]|FList],
        get_prob_handle_nvars(H, NV0),
        ( VList = [E|_], compound(E) ->
            /* multiple tuples of variables */
            ( foreach(Vs, VList), fromto(NV0, NV1,NV2, NV),
              param(H),
              foreach(GVArr, GVs)
            do
                check_collection_to_list(Vs, VsList),
                ec_to_gecode_varlist1(VsList, H, NV1,NV2, GVsList, _),
                GVArr =.. [[]|GVsList]
            )
        ;
            /* single tuple of variables */
            ec_to_gecode_varlist1(VList, H, NV0,NV, GVList, _),
            VArr =.. [[]|GVList],
            GVs = [VArr]
        ), !,
        update_gecode_with_default_newvars(H, NV0, NV),
        g_create_dfa_handle(TArr, Start, FArr, DfaH),
        post_new_event(post_extensional(ConLev, GVs, DfaH), H).
extensional_c(Vars, Transitions, Start, Finals, _ConLev) :-
        get_bip_error(E),
        error(E,extensional(Vars, Transitions, Start, Finals)).


regular(Vars, RegExp) :-
        regular_c(Vars, RegExp, default).

regular_c(Vars, RegExp0, ConLev) :-
        check_compound_to_list(Vars, VList),
        check_regexp(RegExp0, RegExp),
        get_prob_handle_nvars(H, NV0),
        ( VList = [E|_], compound(E) ->
            /* multiple tuples of variables */
            ( foreach(Vs, VList), fromto(NV0, NV1,NV2, NV),
              param(H),
              foreach(GVArr, GVs)
            do
                check_collection_to_list(Vs, VsList),
                ec_to_gecode_varlist1(VsList, H, NV1,NV2, GVsList, _),
                GVArr =.. [[]|GVsList]
            )
        ;
            /* single tuple of variables */
            ec_to_gecode_varlist1(VList, H, NV0,NV, GVList, _),
            VArr =.. [[]|GVList],
            GVs = [VArr]
        ), !,
        update_gecode_with_default_newvars(H, NV0, NV),
        g_create_regdfa_handle(RegExp, DfaH),
        post_new_event(post_extensional(ConLev, GVs, DfaH), H).
regular_c(Vars, RegExp, _ConLev) :-
        get_bip_error(E),
        error(E,regular(Vars, RegExp)).

check_regexp(I0, I) ?-
        integer(I0), !,
        I = I0.
check_regexp(*(E0), E) ?- !,
        E = *(E1),
        check_regexp(E0, E1).
check_regexp(+(E0), E) ?- !,
        E = +(E1),
        check_regexp(E0, E1).
check_regexp(E0+E1, E) ?- !,
        E = E2 + E3,
        check_regexp(E0, E2),
        check_regexp(E1, E3).
check_regexp((E0|E1), E) ?- !,
        E = (E2 | E3),
        check_regexp(E0, E2),
        check_regexp(E1, E3).
check_regexp((E0,{N,M}), E) ?- !, % {N,M} must be before {N} case
        check_integer(N),
        check_integer(M),
        check_nonnegative(N),
        check_nonnegative(M),
        E = (E1,r(N,M)),
        check_regexp(E0, E1).
check_regexp((E0,{N}), E) ?- !,
        check_integer(N),
        check_nonnegative(N),
        E = (E1,{N}),
        check_regexp(E0, E1).
check_regexp(Is, E) :-
        collection_to_list(Is, IList),
        IList = [_|_], !,
        (foreach(I, IList) do check_integer(I)),
        E =.. [[]|IList].
check_regexp(_, _) :-
        set_bip_error(5).

                          
%------------------------------------------------------------------------
% Events and clone management

gfd_update :-
        get_prob_handle(H),
        restore_space_if_needed(H, _SpH),
        H = gfd_prob{nevents:NE, space:Current, last_anc: Anc},
        ( NE == 0 ->
            true % don't update if there are no changes
        ;
            replace_last_ancestor(H, Current, Anc)
        ).

update_space_with_events(H) :-
        H = gfd_prob{events:Es,space:gfd_space{handle:SpH}},
        g_stop_caching(SpH),
        update_space_with_events1(Es, H),
        % no need to check Inst/Chg -- not updated with First=0
        g_propagate_recompute(SpH).


update_space_with_events1(Es, _) :-
        var(Es), !.
%        writeln("*****done updating clone"). 
update_space_with_events1([E|Es], H) ?-
%        write("re"),
        do_event(E, H, _),
        update_space_with_events1(Es, H).
        

% posting an event that may have additional "auxillary" events
post_new_event_with_aux(Es, EsTail, H) :-
        set_new_event(Es, EsTail, H),
        H = gfd_prob{space:gfd_space{handle:SpH}},
        post_new_events1(Es, SpH, DoProp),
        g_start_caching(SpH),
        ( var(DoProp) -> true ; try_propagate(H), wake ).


post_new_events1(Es, _SpH, _DoProp) :-
        var(Es), !.
post_new_events1([E|Es1], SpH, DoProp) ?-
%        writeln(doing-E),
        do_event1(E, SpH, DoProp),
        g_stop_caching(SpH), % first event will be done with caching
        post_new_events1(Es1, SpH, DoProp).


/* post_new_event call wake at the end, so that gfd_do_propagate will
   be woken and executed. This may not be appropriate if the posted 
   event is part of an atomic operation, e.g. in the unify handler, when
   goals only woken later when all the handlers have been executed
*/
post_new_event(E, H) :-
        post_new_event_no_wake(E, H),
        wake.

post_new_event_no_wake(E, H) :-
        Es = [E|ET],
        set_new_event(Es, ET, H),
        do_event(E, H, DoProp), 
        (DoProp == [] -> try_propagate(H) ; true).

try_propagate(H) :-
        H = gfd_prob{prop:Susp},
        schedule_suspensions(1, s([Susp])).


:- demon gfd_do_propagate/1.
gfd_do_propagate(H) :-
        H = gfd_prob{space:gfd_space{handle:SpH},vars:VArr},
        g_propagate(SpH, 1, InstList, ChgList),
        propagate_gecode_changes(SpH, VArr, InstList, ChgList).


propagate_gecode_changes(SpH, VArr, InstList, ChgList) :-
        ( InstList == [] ->
%            store_inc(stats, pg0)
            true
        ;
%            length(InstList, Len), concat_atom([pg, Len], Key), store_inc(stats, Key), 
            ( foreach(Idx, InstList), param(VArr, SpH) do
                arg(Idx, VArr, V),
                (integer(V) -> true ; mark_var_as_set(V),g_get_var_value(SpH, Idx, V))
%                (integer(V) -> true ; g_get_var_value(SpH, Idx, V))
            )
        ),
        ( ChgList == [] ->
            true
        ;
            ( foreach(CIdx, ChgList), param(VArr) do
                arg(CIdx, VArr, U),
                notify_constrained(U),
                get_gecode_attr(U, Attr),
                % assuming only one single problem, otherwise need H
                schedule_suspensions(any of gfd, Attr)
            )
        ).
           

mark_var_as_set(_{gfd{set:S}}) ?- !, S = [].



set_new_event(Es, EsT, H) :-
        check_and_update_handle(H),
        % access H *only* after possible update of handle!
        H = gfd_prob{events_tail:Es,nevents:NE0,space:Sp},
        % Es is joined to events list and events_tail updated to new tail
        setarg(events_tail of gfd_prob, H, EsT), 
        g_trail_undo_for_event(Sp),
        NE1 is NE0+1,
        setarg(nevents of gfd_prob, H, NE1).

% should only be called with a new event
check_and_update_handle(H) :-
        restore_space_if_needed(H, _SpH),
        check_and_update_ancestors(H).


% can be called outside of a new event (e.g. when state is required)
restore_space_if_needed(H, SpH) :-
        H = gfd_prob{space: gfd_space{handle:SpH},last_anc:Anc},
        % pass Anc rather than the C handle, because Anc can be []
        g_check_handle(SpH, Anc, Cloned),
        ( Cloned == [] -> update_space_with_events(H) ; true).


% should only be called with a new event
check_and_update_ancestors(H) :-
        timestamp_age(H, cp_stamp of gfd_prob, Age),
        (Age == old ->
            % first event after a choicepoint
            H = gfd_prob{nevents:NE, events:E,space:Current},
            ( g_state_is_stable(Current) ->
                % only clone if state is stable, i.e. have propagation done
                % this may not be the case if propagation has been delayed
                ( NE >= cloning_distance  ->
                    do_update_ancestor(H, Current)
                 ; E == update ->
                    do_update_ancestor(H, Current)
                ;
                    true
                )
            ;
                true
            ),
            timestamp_update(H, cp_stamp of gfd_prob)
        ;
            H = gfd_prob{nevents:NE, space:Current, last_anc: Anc},
            ( NE > events_max ->
                replace_last_ancestor_if_det(H, Current, Anc)
            ;
                true
            )
        ).

replace_last_ancestor_if_det(H, Current, []) :- !,
        replace_last_ancestor(H, Current, []).
replace_last_ancestor_if_det(H, Current, Anc) :-
         (timestamp_older(Anc, stamp of gfd_space, 
                          Current, stamp of gfd_space) ->
             true
         ;
             replace_last_ancestor(H, Current, Anc)
             % replace ancestor only if computation is deterministic
             % since the ancestor
        ).

% create new ancestor, replacing last ancestor if possible
% this is to reduce amount of subsequent recomputation 
replace_last_ancestor(H, Current, Anc) :-
        ( g_state_is_stable(Current) ->
            do_update_ancestor(H, Current),
            ( Anc = gfd_space{handle:ASpH} ->
                ( timestamp_older(Anc, stamp of gfd_space, 
                                  Current, stamp of gfd_space) ->
                    true 
                ;
                    g_delete(ASpH) % old ancestor replaced and can be deleted
                )
            ;
                true % no old ancestor to delete
            )
        ;
            true
        ).


% clone the current space and make it the last ancestor
do_update_ancestor(H, Current) :-
        new_space_handle(New),
        setarg(last_anc of gfd_prob, H, Current),
        setarg(space of gfd_prob, H, New),
        % setarg/3 used to initialise events list so that Tail is 
        % allocated outside of H 
        setarg(events_tail of gfd_prob, H, Tail),
        setarg(events of gfd_prob, H, Tail),
        setarg(nevents of gfd_prob, H, 0),
        New = gfd_space{handle:NewH},
        g_check_handle(NewH, Current, _Cloned).



/* do_event execute a Gecode event, either for the first time or
   during recomputation. 
   DoProp is set to [] to indicate if propagation should be done after the
   event (used on first execution only).  [] is used as it is slightly
   faster in type indexing.
   Handling of events should not assume that any arguments that are 
   variables in first execution will be variables during recomputation
*/
do_event(E, H, DoProp) :-
%        writeln(doing-E),
        H = gfd_prob{space:gfd_space{handle:SpH}},
        do_event1(E, SpH, DoProp).

do_event1(post_rc(ConLev, GExpr), SpH, DoProp) ?-
        DoProp = [],
        g_post_intrel_cstr(SpH, GExpr, ConLev).
do_event1(post_bool_connectives(ConLev, GBCon), SpH, DoProp) ?-
        DoProp = [],
        g_post_bool_connectives(SpH, GBCon, ConLev).
do_event1(post_alldiff(ConLev, GArray), SpH, DoProp) ?-
        DoProp = [],
        g_post_alldiff(SpH, GArray, ConLev).
do_event1(post_alldiff_offsets(ConLev, GArray, OArray), SpH, DoProp) ?-
        DoProp = [],
        g_post_alldiff_offsets(SpH, GArray, OArray, ConLev).
do_event1(post_nvalues(GArray,Rel,N), SpH, DoProp) ?-
        DoProp = [],
        g_post_nvalues(SpH, GArray, Rel, N).
do_event1(post_count(ConLev, Value,GArray,Rel, N), SpH, DoProp) ?-
        DoProp = [],
        g_post_count(SpH, Value, GArray, Rel, N, ConLev).
do_event1(post_among(Values,GArray,Rel, N), SpH, DoProp) ?-
        DoProp = [],
        g_post_among(SpH, Values, GArray, Rel, N).
do_event1(post_count_matches(Values,GArray,Rel, N), SpH, DoProp) ?-
        DoProp = [],
        g_post_count_matches(SpH, Values, GArray, Rel, N).
do_event1(post_gcc(ConLev, Vals,Occs,GVs), SpH, DoProp) ?-
        DoProp = [],
        g_post_gcc(SpH, Vals, Occs, GVs, ConLev).
do_event1(post_element(ConLev, GI,GArray,GValue), SpH, DoProp) ?-
        DoProp = [],
        g_post_element(SpH, GI, GArray, GValue, ConLev).
do_event1(post_sequence(ConLev, Lo, Hi, K, VarArray, ValArray), SpH, DoProp) ?-
        DoProp = [],
        g_post_sequence(SpH, Lo, Hi, K, VarArray, ValArray, ConLev).
do_event1(post_sequence_01(ConLev, Lo, Hi, K, VarArray), SpH, DoProp) ?-
        DoProp =[],
        g_post_sequence_01(SpH, Lo, Hi, K, VarArray, ConLev).
do_event1(post_sorted2(ConLev, UsArray, SsArray), SpH, DoProp) ?-
        DoProp = [],
        g_post_sorted2(SpH, UsArray, SsArray, ConLev).
do_event1(post_sorted(ConLev, UsArray, SsArray, PsArray), SpH, DoProp) ?-
        DoProp = [],
        g_post_sorted(SpH, UsArray, SsArray, PsArray, ConLev).
do_event1(post_circuit(ConLev, SArray, Offset), SpH, DoProp) ?-
        DoProp = [],
        g_post_circuit(SpH, SArray, Offset, ConLev).
do_event1(post_circuit_cost(ConLev, SArray, CMArray, ACArray, GCost, Offset), SpH, DoProp) ?-
        DoProp = [],
        g_post_circuit_cost(SpH, SArray, CMArray, ACArray, GCost, Offset, ConLev).
do_event1(post_ham_path(ConLev,Start,End,SArray, Offset), SpH, DoProp) ?-
        DoProp = [],
        g_post_ham_path(SpH, Start, End, SArray, Offset, ConLev).
do_event1(post_ham_path_cost(ConLev,Start,End,SArray,CMArray,ACArray,GCost,Offset), SpH, DoProp) ?-
        DoProp = [],
        g_post_ham_path_cost(SpH, Start, End, SArray, CMArray, ACArray, GCost, Offset, ConLev).
do_event1(post_disj(StartArray,DurArray,SchArray), SpH, DoProp) ?-
        DoProp = [],
        % ConLev not supported for this constraint
        g_post_disj(SpH, StartArray, DurArray, SchArray).
do_event1(post_disjflex(StartArray,DurArray,EndArray,SchArray), SpH, DoProp) ?-
        DoProp = [],
        % ConLev not supported for this constraint
        g_post_disjflex(SpH, StartArray, DurArray, EndArray, SchArray).
do_event1(post_cumulatives(_ConLev, Starts,Durations,Ends,Usages,Used,Limits,AtMost), SpH, DoProp) ?-
        % ignore ConLev for now (only gfd_vc allowed)
        DoProp = [],
        g_post_cumulatives(SpH, Starts, Durations, Ends, Usages, Used, Limits, AtMost).
do_event1(post_cumulative(Starts,Durations,Usages,Limit,Schs), SpH, DoProp) ?-
        DoProp = [],
        g_post_cumulative(SpH, Starts, Durations, Usages, Limit, Schs).
do_event1(post_cumulativeflex(Starts,Durations,Ends,Usages,Limit,Schs), SpH, DoProp) ?-
        DoProp = [],
        g_post_cumulativeflex(SpH, Starts, Durations, Ends, Usages, Limit, Schs).
do_event1(post_disjoint2(Xs,Widths,Ys,Heights,Optionals), SpH, DoProp) ?-
        DoProp = [],
        g_post_disjoint2(SpH, Xs, Widths, Ys, Heights, Optionals).
do_event1(post_disjointflex2(X1s,Widths,Y1s,Heights,Optionals,X2s,Y2s), SpH, DoProp) ?-
        DoProp = [],
        g_post_disjointflex2(SpH, X1s, Widths, Y1s, Heights, Optionals, X2s, Y2s).
do_event1(post_precede(S, T, GArray), SpH, DoProp) ?-
        DoProp = [],
        g_post_precede(SpH, S, T, GArray).
do_event1(post_precede_chain(Vals, GArray), SpH, DoProp) ?-
        DoProp = [],
        g_post_precede_chain(SpH, Vals, GArray).
do_event1(post_interval(GArray,Lo,Hi), SpH, DoProp) ?-
        DoProp = [],
        g_post_interval(SpH, GArray, Lo, Hi).
do_event1(post_var_interval_reif(GV,Lo,Hi,GBool), SpH, DoProp) ?-
        DoProp = [],
        g_post_var_interval_reif(SpH, GV, Lo, Hi, GBool).
do_event1(post_dom(GArray,DArray), SpH, DoProp) ?-
        DoProp = [],
        g_post_dom(SpH, GArray, DArray).
do_event1(post_dom_handle(GArray,DArray), SpH, DoProp) ?-
        DoProp = [],
        g_post_dom_handle(SpH, GArray, DArray).
do_event1(post_var_dom_reif(GV,DArray,GBool), SpH, DoProp) ?-
        DoProp = [],
        g_post_var_dom_reif(SpH, GV, DArray, GBool).
do_event1(post_var_val_reif(GV,Value,GBool), SpH, DoProp) ?-
        DoProp = [],
        g_post_var_val_reif(SpH, GV, Value, GBool).
do_event1(post_exclude_dom(GArray,DArray), SpH, DoProp) ?-
        DoProp = [],
        g_post_exclude_dom(SpH, GArray, DArray).
do_event1(post_exclude_dom_handle(GArray,Dom), SpH, DoProp) ?-
        DoProp = [],
        g_post_exclude_dom_handle(SpH, GArray, Dom).
do_event1(post_exclude_val(GArray,Val), SpH, DoProp) ?-
        DoProp = [],
        g_post_exclude_val(SpH, GArray, Val).
do_event1(post_exclude_range(GArray,Lo,Hi), SpH, DoProp) ?-
        DoProp = [],
        g_post_exclude_range(SpH, GArray, Lo, Hi).
do_event1(newvars_interval(NV,Lo,Hi), SpH, _DoProp) ?-
        %DoProp = 0,
        g_add_newvars_interval(SpH, NV, Lo, Hi).
do_event1(newvars_dom(NV,DArray), SpH, _DoProp) ?-
        %DoProp = 0,
        g_add_newvars_dom(SpH, NV, DArray).
do_event1(newvars_dom_handle(NV,Dom), SpH, _DoProp) ?-
        %DoProp = 0,
        g_add_newvars_dom_handle(SpH, NV, Dom).
do_event1(newvars_dom_union(GX,GY,NV), SpH, _DoProp) ?-
        %DoProp = 0,
        g_add_newvars_dom_union(SpH, NV, GX, GY).
do_event1(newboolvars(NV,VArr), SpH, _DoProp) ?-
        %DoProp = 0,
        g_add_newvars_as_bool(SpH, NV, VArr).
do_event1(connectnewbools(VArr), SpH, _DoProp) ?-
        %DoProp = 0,
        g_link_newbools(SpH, VArr).
do_event1(copyvar(NV,OldIdx), SpH, _DoProp) ?-
        %DoProp = 0,
        g_add_newvar_copy(SpH, NV, OldIdx).
do_event1(setvar(Idx, Val), SpH, DoProp) ?-
        DoProp = [],
        g_post_setvar(SpH, Idx, Val).
do_event1(post_exclude_var_val(Idx, Val), SpH, DoProp) ?-
        DoProp = [],
        g_post_exclude_var_val(SpH, Idx, Val).
do_event1(post_sum(ConLev, GArray, Rel, C), SpH, DoProp) ?-
        DoProp = [],
        g_post_sum(SpH, GArray, Rel, C, ConLev).
do_event1(post_sum_reif(ConLev, GArray, Rel, C, GBool), SpH, DoProp) ?-
        DoProp = [],
        g_post_sum_reif(SpH, GArray, Rel, C, GBool, ConLev).
do_event1(post_lin(ConLev, GArray, CArray, Rel, C), SpH, DoProp) ?-
        DoProp = [],
        g_post_lin(SpH, GArray, CArray, Rel, C, ConLev).
do_event1(post_lin_reif(ConLev, GArray, CArray, Rel, C, GBool), SpH, DoProp) ?-
        DoProp = [],
        g_post_lin_reif(SpH, GArray, CArray, Rel, C, GBool, ConLev).
do_event1(post_maxlist(ConLev, GV, GArray), SpH, DoProp) ?-
        DoProp = [],
        g_post_maxlist(SpH, GV, GArray, ConLev).
do_event1(post_minlist(ConLev, GV, GArray), SpH, DoProp) ?-
        DoProp = [],
        g_post_minlist(SpH, GV, GArray, ConLev).
do_event1(post_sqrt(ConLev, GRes, GX), SpH, DoProp) ?-
        DoProp = [],
        g_post_sqrt(SpH, GRes, GX, ConLev).
do_event1(post_sq(ConLev, GRes, GX,_Power), SpH, DoProp) ?-
        DoProp = [],
        g_post_sq(SpH, GRes, GX, ConLev).
do_event1(post_sq(ConLev, GRes, GX), SpH, DoProp) ?-
        DoProp = [],
        g_post_sq(SpH, GRes, GX, ConLev).
do_event1(post_abs(ConLev, GRes, GX), SpH, DoProp) ?-
        DoProp = [],
        g_post_abs(SpH, GRes, GX, ConLev).
do_event1(post_div(_ConLev, GRes, GX, GY), SpH, DoProp) ?-
        DoProp = [],
        g_post_div(SpH, GRes, GX, GY).
do_event1(post_mult(ConLev, GRes, GX, GY), SpH, DoProp) ?-
        DoProp = [],
        g_post_mult(SpH, GRes, GX, GY, ConLev).
do_event1(post_mod(_ConLev, GRes, GX, GY), SpH, DoProp) ?-
        DoProp = [],
        g_post_mod(SpH, GRes, GX, GY).
do_event1(post_divmod(_ConLev, GX, GY, GQ, GM), SpH, DoProp) ?-
        DoProp = [],
        g_post_divmod(SpH, GX, GY, GQ, GM).
do_event1(post_min2(ConLev, GRes, GX, GY), SpH, DoProp) ?-
        DoProp = [],
        g_post_min2(SpH, GRes, GX, GY, ConLev).
do_event1(post_max2(ConLev, GRes, GX, GY), SpH, DoProp) ?-
        DoProp = [],
        g_post_max2(SpH, GRes, GX, GY, ConLev).
do_event1(post_rel(ConLev, GX, Rel, GY), SpH, DoProp) ?-
        DoProp = [],
        g_post_rel(SpH, GX, Rel, GY, ConLev).
do_event1(post_collection_rel(ConLev, GXs, Rel, GY), SpH, DoProp) ?-
        DoProp = [],
        g_post_collection_rel(SpH, GXs, Rel, GY, ConLev).
do_event1(post_lwb(GV, Lwb), SpH, DoProp) ?-
        DoProp = [],
        g_post_lwb(SpH, GV, Lwb).
do_event1(post_upb(GV, Lwb), SpH, DoProp) ?-
        DoProp = [],
        g_post_upb(SpH, GV, Lwb).
do_event1(newbool(Idx,BIdx), SpH, _DoProp) ?-
        %DoProp = 0,
        g_add_newbool(SpH, Idx, BIdx).
do_event1(post_boolchannel(ConLev,GV,GBArr,Min), SpH, DoProp) ?-
        DoProp = [],
        g_post_boolchannel(SpH, GV, GBArr, Min, ConLev).
do_event1(post_inverse(ConLev,Arr1,Arr2), SpH, DoProp) ?-
        DoProp = [],
        g_post_inverse(SpH, Arr1, Arr2, ConLev).
do_event1(post_inverse_offset(ConLev,Arr1,Off1,Arr2,Off2), SpH, DoProp) ?-
        DoProp = [],
        g_post_inverse_offset(SpH, Arr1, Off1, Arr2, Off2, ConLev).
do_event1(post_ordered(ConLev,XArr,Rel), SpH, DoProp) ?-
        DoProp = [],
        g_post_ordered(SpH, XArr, Rel, ConLev).
do_event1(post_lex_order(ConLev,XArr,Rel,YArr), SpH, DoProp) ?-
        DoProp = [],
        g_post_lex_order(SpH, XArr, Rel, YArr, ConLev).
do_event1(post_bin_packing(IArr,SArr,LArr), SpH, DoProp) ?-
        DoProp = [],
        g_post_bin_packing(SpH, IArr, SArr, LArr).
do_event1(post_mem(VArr,Mem), SpH, DoProp) ?-
        DoProp = [],
        g_post_mem(SpH, VArr, Mem).
do_event1(post_mem_reif(VArr,Mem, Bool), SpH, DoProp) ?-
        DoProp = [],
        g_post_mem_reif(SpH, VArr, Mem, Bool).
do_event1(post_table(ConLev,VArr,TsH,Size,Emph), SpH, DoProp) ?-
        DoProp = [],
        g_post_table(SpH, VArr, TsH, Size, Emph, ConLev).
do_event1(post_extensional(ConLev,GVs,DfaH), SpH, DoProp) ?-
        DoProp = [],
        g_post_extensional(SpH, GVs, DfaH, ConLev).



%------------------------------------------------------------------------
% labelling and search

indomain(V) :- do_indomain_min(V).


indomain(V, min) ?- !,
        do_indomain_min(V).
indomain(V, max) ?- !,
        do_indomain_max(V).
indomain(V, median) ?- !,
        do_indomain_median(V).
indomain(V, middle) ?- !,
        do_indomain_middle(V).
indomain(V, enum) ?- !,
        indomain(V).
indomain(V, I) ?- 
        integer(I), !,
        do_indomain_from_value(V, I).
indomain(V, I) ?-
        error(6, indomain(V,I)).


labeling(Vars) :-
        check_collection_to_list(Vars, List), !,
        gfd_update,
        ( foreach(Var, List) do do_indomain_min(Var) ).
labeling(Vars) :-
        error(5, labeling(Vars)).

labeling(Vars, Select, Choice) :-
        select_var_setup(Vars, 0, H, VsH), !,
        gfd_update,
        labeling1(VsH, H, Select, Choice).
labeling(Vars, Select, Choice) :-
        error(5, labeling(Vars, Select, Choice)).


labeling1(VsH, H, Select, Choice) :-
        H = gfd_prob{vars:Vars,space:gfd_space{handle:SpH}},
        (g_select(SpH, VsH, Select, Idx) ->
            arg(Idx, Vars, V),
            try_value1(V, SpH, H, Idx, Choice),
            labeling1(VsH, H, Select, Choice)
        ;
            true
        ).

                             
do_indomain_min(I) :- integer(I), !.
do_indomain_min(V{gfd:Attr}) ?-
        nonvar(Attr), !,
        Attr = gfd{prob:H, idx:Idx},
        restore_space_if_needed(H, SpH),
        g_get_var_lwb(SpH, Idx, Lo),
        % -1 for lower bound
        indomain_and_prop(Idx, V, H, -1, Lo).
do_indomain_min(V) :-
        error(5, indomain(V, min)).


do_indomain_max(I) :- integer(I), !.
do_indomain_max(V{gfd:Attr}) ?-
        nonvar(Attr), !,
        Attr = gfd{prob:H, idx:Idx},
        restore_space_if_needed(H, SpH),
        g_get_var_upb(SpH, Idx, Lo),
        % 1 for upper bound
        indomain_and_prop(Idx, V, H, 1, Lo).
do_indomain_max(V) :-
        error(5, indomain(V, max)).


% this does not create a useless choice-point for the last element
% in the domain. This avoids the wasteful recreation of a valid space
% for the backtrack beyond the last element, only to fail immediately.
indomain_and_prop(_Idx, V, _H, _W, V) ?-
        integer(V). % last element of domain
indomain_and_prop(Idx, V, H, Which, Val) :- 
        meta(V),
        ( V = Val
        ;
          restore_space_if_needed(H, _SpH), 
          % -1 for lower bound
          arg(space of gfd_prob, H, Sp),
          g_update_and_get_var_bound(Sp, Idx, Val, Which, Bound),
          indomain_and_prop(Idx, V, H, Which, Bound)
        ).


% The Gecode-style value choice, not necessarily grounding the variable
% The indomain style value choice remove old values directly (i.e. not
% via events), and so is not performed during recomputation. This
% means that there is one event per variable in the search, so the
% total number of events are bounded by the number of variables. 
try_value(X, _Method) :- integer(X).
try_value(X{gfd:Attr}, Method) ?- 
        nonvar(Attr),
        Attr = gfd{idx:Idx,prob:H},
        restore_space_if_needed(H, SpH),
	try_value1(X, SpH, H, Idx, Method).

try_value1(X, SpH, H, Idx, min) :-
        g_get_var_lwb(SpH, Idx, Val),
        ( X=Val 
        ;                 
          post_new_event(post_exclude_var_val(Idx, Val), H)
        ).
try_value1(X, SpH, H, Idx, max) :-
        g_get_var_upb(SpH, Idx, Val),
        ( X=Val 
        ;                 
          post_new_event(post_exclude_var_val(Idx, Val), H)
        ).
try_value1(X, SpH, H, Idx, median) :-
        g_get_var_median(SpH, Idx, Val),
        ( X=Val 
        ;                 
          post_new_event(post_exclude_var_val(Idx, Val), H)
        ).
try_value1(_X, SpH, H, Idx, split) :-
	g_get_var_bounds(SpH, Idx, Lo, Hi),
	Split is (Lo+Hi) div 2,
        gfdvar(Idx, _, GX),
        ( post_new_event(post_rel(gfd_bc, GX, (#=<), Split), H)
        ;
          post_new_event(post_rel(gfd_bc, GX, (#>), Split), H)
        ).
try_value1(_X, SpH, H, Idx, reverse_split) :-
	g_get_var_bounds(SpH, Idx, Lo, Hi),
	Split is (Lo+Hi) div 2,
        gfdvar(Idx, _, GX),
        ( post_new_event(post_rel(gfd_bc, GX, (#>), Split), H)
        ;
          post_new_event(post_rel(gfd_bc, GX, (#=<), Split), H)
        ).
% indomain style value-choice
/* version that does remove last value on backtracking
try_value1(X, SpH, H, Idx,indomain_min) :-
        g_get_var_lwb(SpH, Idx, Lo),
        ( X = Lo
        ;
          post_new_event(post_exclude_var_val(Idx, Lo), H),
          H = gfd_prob{space:gfd_space{handle:SpH1}},
          (integer(X) -> 
              true
          ;
              try_value1(X, SpH1, H, Idx,indomain_enum)
          )
        ).
*/
try_value1(X, SpH, H, Idx, indomain_min) :-
        g_get_var_lwb(SpH, Idx, Lo),
        % -1 for lower bound
        indomain_and_prop(Idx, X, H, -1, Lo).
try_value1(X, SpH, H, Idx, indomain_max) :-
        g_get_var_upb(SpH, Idx, Hi),
        % -1 for lower bound
        indomain_and_prop(Idx, X, H, 1, Hi).
try_value1(X, SpH, H, Idx, indomain_median) :-
        g_get_var_median(SpH, Idx, Med),
        arg(space of gfd_prob, H, Sp),
        indomain_from(X, Med, H, Sp, Idx).
try_value1(X, SpH, H, Idx, indomain_middle) :-
        g_get_var_bounds(SpH, Idx, Lo, Hi),
        Mid is (Hi+Lo) // 2,
        arg(space of gfd_prob, H, Sp),
        indomain_from(X, Mid, H, Sp, Idx).
try_value1(X, _SpH, H, Idx, indomain_from(I)) :-
        integer(I),
        arg(space of gfd_prob, H, Sp),
        indomain_from(X, I, H, Sp, Idx).

% Possible extensions
%try_value1(X, middle) :-
%try_value1(X, random) :-
%try_value1(X, indomain_interval) :-
%try_value1(X, indomain_interval_min) :-
%try_value1(X, indomain_interval_max) :-


do_indomain_median(I) :- integer(I), !.
do_indomain_median(V{gfd:Attr}) ?-
        nonvar(Attr), !,
        Attr = gfd{prob:H, idx:Idx},
        restore_space_if_needed(H, SpH),
        g_get_var_median(SpH, Idx, Med),
        arg(space of gfd_prob, H, Sp),
        indomain_from(V, Med, H, Sp, Idx).
do_indomain_median(V) :-
        error(5, indomain(V, median)).


do_indomain_middle(I) :- integer(I), !.
do_indomain_middle(V{gfd:Attr}) ?-
        nonvar(Attr), !,
        Attr = gfd{prob:H, idx:Idx},
        restore_space_if_needed(H, SpH),
        g_get_var_bounds(SpH, Idx, Lo, Hi),
        Mid is (Hi+Lo) // 2,
        arg(space of gfd_prob, H, Sp),
        indomain_from(V, Mid, H, Sp, Idx).
do_indomain_middle(V) :-
        error(5, indomain(V, middle)).


do_indomain_from_value(I, Value) :- 
        integer(I), !,
        Value == I.
do_indomain_from_value(V{gfd:Attr}, Value) ?-
        nonvar(Attr), !,
        Attr = gfd{prob:H, idx:Idx},
        restore_space_if_needed(H, _SpH),
        arg(space of gfd_prob, H, Sp),
        indomain_from(V, Value, H, Sp, Idx).
do_indomain_from_value(V,I) :-
        error(5, indomain(V, I)).


/* label V starting with Val, and then trying values next to Val,
   alternating between higher and lower values. OldHi and OldLo
   are used to record the last Hi and Lo values tried.
*/
indomain_from(V, Val, H, Sp, Idx) :- 
        Val0 is Val + 1,
        shelf_create(last(Val), OldHi),
        shelf_create(last(Val0), OldLo),
        indomain_from1(1, V, H, Sp, Idx, OldHi, OldLo).

indomain_from1(-1, V, H, Sp, Idx, OldHi, OldLo) :-
        % trying larger values
        shelf_get(OldHi, 1, Old),
        ( g_update_and_get_var_bound(Sp, Idx, Old, -1, Next) ->
            % propagate OldHi value as lower bound, to get next lowest
            % value in domain to try
            shelf_set(OldHi, 1, Next),
            V = Next
        ;
            !,
            % no higher values left in domain, try smaller values only
            shelf_get(OldLo, 1, Lo0),
            restore_space_if_needed(H, _SpH),
            g_update_and_get_var_bound(Sp, Idx, Lo0, 1, Lo), % may fail!
            indomain_and_prop(Idx, V, H, 1, Lo)
        ).
indomain_from1(1, V, H, Sp, Idx, OldHi, OldLo) :-
        % trying smaller values
        shelf_get(OldLo, 1, Old),
        ( g_update_and_get_var_bound(Sp, Idx, Old, 1, Next) ->
            shelf_set(OldLo, 1, Next),
            V = Next
        ;
            !,
            shelf_get(OldHi, 1, Hi0),
            restore_space_if_needed(H, _SpH),
            g_update_and_get_var_bound(Sp, Idx, Hi0, -1, Hi), % may fail
            indomain_and_prop(Idx, V, H, -1, Hi)
        ).
indomain_from1(Which, V, H, Sp, Idx, OldHi, OldLo) :-
        restore_space_if_needed(H, _SpH),
        NewWhich is Which * -1,
        indomain_from1(NewWhich, V, H, Sp, Idx, OldHi, OldLo).


select_var(X, Xs, IdxsH, Arg, Select) :-
        ( is_handle(Xs) ->
            Xs = IdxsH,
            get_prob_handle(H),
	    select_var1(X, H, Select, IdxsH)
        ;
            select_var_setup(Xs, Arg, H, IdxsH),
	    select_var1(X, H, Select, IdxsH)
        ).


select_var_setup(Xs, Arg, H, IdxsH) :-
        check_collection_to_list(Xs, LXs),
        ( foreach(V0, LXs), 
          param(Arg),
          fromto(LIdxs, LIdxs0,LIdxs1, [])
        do
            (Arg == 0 -> V = V0 ; arg(Arg, V, V0)), 
            (get_gecode_var(V, GV), gfdvar(VIdx,_,GV) ->
                LIdxs0 = [VIdx|LIdxs1]
            ;
                LIdxs0 = LIdxs1
            )
        ),
        get_prob_handle(H),
        Idxs =.. [[]|LIdxs],
        g_create_idxs_handle(Idxs, IdxsH).

select_var1(X, H, Select, IdxsH) :-
        restore_space_if_needed(H, SpH),
        g_select(SpH, IdxsH, Select, XIdx),  % May fail
        H = gfd_prob{vars:PVs},
        arg(XIdx, PVs, X).


% Selection criteria evaluators for gfd_search's delete/5 and search/6

:- export
	max_regret_lwb/2,
	max_regret_upb/2,
	max_weighted_degree/2,
	most_constrained_per_value/2,
	max_weighted_degree_per_value/2.

max_regret_lwb(X, Number) :-
	( nonvar(X) -> true ; Number is -get_regret_lwb(X) ).

max_regret_upb(X, Number) :-
	( nonvar(X) -> true ; Number is -get_regret_upb(X) ).

max_weighted_degree(X, Number) :-
	( nonvar(X) -> true ; Number is -get_weighted_degree(X) ).

most_constrained_per_value(X, Number) :-
	( nonvar(X) ->
	    true	% pick constants first and commit
	;
            get_constraints_number(X, NC),
            get_domain_size(X, Size),
            Number is fix(round(Size/NC))
        ).

max_weighted_degree_per_value(X, Number) :-
	( nonvar(X) ->
	    true	% pick constants first and commit
	;
            get_weighted_degree(X, AFC),
            get_domain_size(X, Size),
            Number is fix(round(Size/AFC))
        ).

                                  
/* search/6 maps to Gecode's branching and search engines. 
   As cloning is done before and after the use of the search engine,
   we cannot treat the search as a normal gfd event. Instead, the 
   event queue is marked, so that the space returned after the search
   will become an ancestor if further events are posted. 
*/   
:- export struct(gfd_stats(prop,fail,nodes,depth,mem)).
:- export struct(gfd_control(commit_distance,adaptive_distance,threads)).

search(Vars, Pos, Select, Choice, Method, Option) :-
        check_collection_to_list(Vars, List),
        check_atom(Select), 
        (var_selection(Select) -> true ; set_bip_error(6)),
        check_atom(Choice), 
        (val_choice(Choice) -> true ; set_bip_error(6)),
        check_integer(Pos),
        ( Pos > 0 ->
            ( foreach(E, List), foreach(V, VList),
              param(Pos)
            do
                arg(Pos, E, V)
            )
        ;
            VList = List
        ),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_varlist1(VList, H, NV0,NV, GList, _),
        (translate_search_method(Method, TMethod) -> true ; set_bip_error(6)),
        GArray =.. [[]|GList],
        process_options(Option, TieBreak, Stats, Stop, Timeout, Control),
        !,
        update_gecode_with_default_newvars(H, NV0, NV),
        restore_space_if_needed(H, SpH),
        H = gfd_prob{space:SP0},
        g_setup_search(SP0, GArray, Select, Choice, TMethod, TieBreak, 
                       Timeout, Stop, Control, EngH),
        new_space_handle(SP),
        setarg(space of gfd_prob, H, SP),
        do_search(SpH, EngH, H, Method, Stats).
search(Vars, Pos, Select, Choice, Method, Option) :-
        get_bip_error(E),
        error(E, search(Vars, Pos, Select, Choice, Method, Option)).


process_options([], TieBreak, _Stats, Stop, Timeout, Control) ?- !,
        (var(TieBreak) -> TieBreak = none ; true),
        (var(Stop) ->  Stop = none ; true),
        (var(Timeout) -> Timeout = 0 ; true),
        (var(Control) -> Control = none ; true).
process_options([O|Os], TieBreak, Stats, Stop, Timeout, Control) :-
        process_option(O, TieBreak, Stats, Stop, Timeout, Control),
        process_options(Os, TieBreak, Stats, Stop, Timeout, Control).

process_option(tiebreak(TieBreak0), TieBreak, _, _, _, _) ?-
        atom(TieBreak0), !,
        var_selection(TieBreak0),
        TieBreak0 = TieBreak.
process_option(stats(Stats0), _, Stats, _, _, _) ?- !,
        (nonvar(Stats0) ->  Stats0 = gfd_stats{} ; true),
        Stats = Stats0.
process_option(backtrack(BT), _, Stats, _, _, _) ?- 
        var(BT), !,
        Stats = gfd_stats{fail:BT}.
process_option(limits(Stop0), _, _, Stop, _, _) ?-
        nonvar(Stop0),
        Stop0 = gfd_stats{}, !,
        Stop = Stop0.
process_option(nodes(N), _, _, Stop, _, _) ?- % ic/fd/gfd_search compatibility
        integer(N), N >= 0, !,
        Stop = gfd_stats{nodes:N}.
process_option(timeout(T), _, _, _, Timeout, _) ?-
        number(T), T >= 0, !,
        convert_timeout_to_ms(T, Timeout).
process_option(control(Control0), _, _, _, _, Control) ?-
        nonvar(Control0), 
        Control0 = gfd_control{}, !,
        Control0 = Control.
process_option(O, _, _, _, _, _) :-
        printf(error, "Unsupported search option %w%n", [O]),
        set_bip_error(6).


convert_timeout_to_ms(T0, T) :-
        T is integer(round(T0*1000)).

:- mode var_selection(+).
var_selection(input_order).
var_selection(first_fail).
var_selection(anti_first_fail).
var_selection(occurrence).
var_selection(anti_occurrence).
var_selection(largest).
var_selection(smallest).
var_selection(largest_uwb).
var_selection(smallest_upb).
var_selection(most_constrained).
var_selection(most_constrained_per_value).
var_selection(least_constrained_per_value).
var_selection(max_regret).
var_selection(max_regret_lwb).
var_selection(min_regret_lwb).
var_selection(max_regret_upb).
var_selection(min_regret_upb).
var_selection(random).
% newer heuristics, added since gecode 3.1
var_selection(max_weighted_degree).
var_selection(min_weighted_degree).
var_selection(max_weighted_degree_per_value).
var_selection(min_weighted_degree_per_value).


:- mode val_choice(+).
val_choice(indomain).
val_choice(indomain_reverse_enum).
val_choice(indomain_min).
val_choice(indomain_max).
val_choice(indomain_median).
val_choice(indomain_random).
val_choice(indomain_split).
val_choice(indomain_reverse_split).
val_choice(indomain_interval).
val_choice(indomain_interval_min).
val_choice(indomain_interval_max).

:- mode translate_search_method(+, -).
translate_search_method(complete, complete).
translate_search_method(lds(D), lds(D)) :- integer(D).
translate_search_method(bb_min(CV),  bb_min(CIdx)) :-
        get_gecode_attr(CV, gfd{idx:CIdx}).
translate_search_method(restart_min(CV),  restart_min(CIdx)) :-
        get_gecode_attr(CV, gfd{idx:CIdx}).

:- mode optimising_search_method(+).
optimising_search_method(bb_min(_)).
optimising_search_method(restart_min(_)).

mark_handle_for_ancestor_update(H) :-
        setarg(events of gfd_prob, H, update).

do_search1(SP, EngH, LastSpH, InstList, ChgList, Stats) :-
        g_do_search(SP, EngH, LastSpH, InstList, ChgList, Stats, Status),
        check_search_status(Status).

check_search_status(0) ?- !, fail.  % search succeeded
check_search_status(1) ?- !.        % search failed
check_search_status(2) ?- !, abort. % search aborted with problem
check_search_status(Status) :-
        Status > 2,                 % search aborted because limit(s) reached 
        report_search_limits(Status),
        printf(log_output, "Search not completed.%n", []),
        (Status /\ 1 =:= 1 -> true/*has sol*/ ; fail/*no sol*/).

% Status here must correspond to those defined in Cutoff's enum in gfd.hpp
report_search_limits(Status) :- 
        ( Status /\  4 =:= 4 ->  % 1<<2
            printf(log_output, "Node limit reached. ", [])
        ;
            true
        ),
        ( Status /\  8 =:= 8  ->  % 1<<3
            printf(log_output, "Failure limit reached. ", [])
        ;
            true
        ),
        ( Status /\ 16 =:= 16 ->  % 1<<4
            printf(log_output, "Time limit reached. ", [])
        ;
            true
        ),
        ( Status /\ 32 =:= 32 ->  % 1<<5
            printf(log_output, "Memory limit reached. ", [])
        ;
            true
        ).
        
do_search(LastSpH, EngH, H, Method, Stats) :-
        mark_handle_for_ancestor_update(H),
        H = gfd_prob{space:SP},
        repeat,
        (do_search1(SP, EngH, LastSpH, InstList, ChgList, Stats) -> % fails if no more solution
            % a clone is created automatically by gecode's search
            % make sure an ancestor will be created on new event
%            writeln(succ:g_do_search(SP, EngH, MethodCode, LastSpH, InstList)),
            % no other solution if optimising
            ( optimising_search_method(Method) -> !  ; true), 
            H = gfd_prob{vars:VArr},
            SP = gfd_space{handle:SpH},
            propagate_gecode_changes(SpH, VArr, InstList, ChgList),
            wake % need wake here -- search is not an event
        ;
%            writeln(failed:g_do_search(SP, EngH, MethodCode, LastSpH, InstList)),
            !, fail
        ).


%------------------------------------------------------------------------
% meta support


% get_min, get_max, get_bounds, get_median -- these are designed for
% value choice for a variable. The following behaviour is implemented:
%  a) If "variable" is an integer, then return its value
%  b) If variable is a (gfd) domain variable, the appropriate
%     value is obtained from Gecode and returned
%  c) If variable is a non-domain variable, turn it into a domain
%     variable with default bounds, and return the appropriate value
%  d) Otherwise, raise a type error
%
% Based on the semantics of equivalent ic predicates

get_min(I, Lo) :- 
        integer(I), !,
        Lo = I.
get_min(V, Lo) :-
        var(V), !,
        get_prob_handle(H),
        ( get_gecode_attr(V, Attr) ->
            Attr = gfd{idx:Idx},
            restore_space_if_needed(H, SpH),
            g_get_var_lwb(SpH, Idx, Lo)
        ;
            create_and_add_default_gfdvar(V, H),
            gfd_get_default(interval_min, Lo)
        ).
get_min(V, Lo) :-
        error(5, get_min(V, Lo)).


get_max(I, Hi) :- 
        integer(I), !,
        Hi = I.
get_max(V, Hi) :-
        var(V), !,
        get_prob_handle(H),
        ( get_gecode_attr(V, Attr) ->
            Attr = gfd{idx:Idx},
            restore_space_if_needed(H, SpH),
            g_get_var_upb(SpH, Idx, Hi)
        ;
            create_and_add_default_gfdvar(V, H),
            gfd_get_default(interval_max, Hi)
        ).
get_max(V, Lo) :-
        error(5, get_max(V, Lo)).


get_bounds(I, Lo, Hi) :-
        integer(I), !,
        Lo = I,
        Hi = I.
get_bounds(V, Lo, Hi) :-
        var(V), !, 
        ( gfd_get_var_bounds(V, Lo, Hi) ->
            true
        ;
            % not an existing gfd var, add it
            get_prob_handle(H),
            create_and_add_default_gfdvar(V, H),
            gfd_default_interval(Lo, Hi)
        ).
get_bounds(V, Lo, Hi) :-
        error(5, get_bounds(V, Lo, Hi)).


get_integer_bounds(V, Lo, Hi) :- % ic compatibility
        get_bounds(V, Lo, Hi).

get_finite_integer_bounds(V, Lo, Hi) :- % ic compatibility
        get_bounds(V, Lo, Hi).


get_median(I, M) :-
        integer(I), !,
        M = I.
get_median(V, Med) :-
        var(V), !,
        get_prob_handle(H),
        ( get_gecode_attr(V, Attr) ->
            Attr = gfd{idx:Idx},
            restore_space_if_needed(H, SpH),
            g_get_var_median(SpH, Idx, Med)
        ;
            H = gfd_prob{nvars:NV0},
            new_gfdvar(V, H, NV0,NV, GV),
            gfdvar(Idx,_BI,GV),
            gfd_default_interval(Min, Max),
            do_update_newvars_with_domain_interval(H, NV, Min, Max),
            % must have a valid space after adding new var
            H = gfd_prob{space:gfd_space{handle:SpH}},
            g_get_var_median(SpH, Idx, Med)
        ).
get_median(V, M) :-
        error(5, get_median(V, M)).

% get_domain/2 does not create a gfd domain variable, instead a type
% error is raised. This is done to be compatible with ic
get_domain(I, Dom) :-
        integer(I), !,
        Dom = [I].
get_domain(_{gfd:Attr}, Dom) ?-
        nonvar(Attr), 
        Attr = gfd{prob:H, idx:Idx}, !,
        restore_space_if_needed(H, SpH),
        g_get_var_domain(SpH, Idx, Dom).
get_domain(X, Dom) :-
        error(5, get_domain(X, Dom)).


get_domain_as_list(V, DomList) :-
        get_domain(V, Dom),
        translate_domain_to_list(Dom, DomList).

  translate_domain_to_list([], DomList) ?- !,
        DomList = [].
  translate_domain_to_list([X|Xs], DomList) ?-
        ( integer(X) ->
            DomList = [X|DomList0]
        ; X = Lo..Hi ->
            ( for(I, Lo, Hi), 
              fromto(DomList, [I|DomList1],DomList1, DomList0)
            do
                true
            )
        ;
            fail
        ),
        translate_domain_to_list(Xs, DomList0).


% get_domain_size, get_delta, get_constraints_number, get_weighted_degree
% get_regret_lwb, get_regret_upb
% These predicates are designed for use with variable selection. They
% follow the common behaviour:
%  a) If "variable" is an integer, return an appropriate value
%  b) If variable is a (gfd) domain variable, obtain the appropriate
%     value from Gecode for this variable
%  c) If variable is not a gfd domain variable, return an appropriate
%     value, but do *not* turn it into a gfd domain variable
%  d) Otherwise, raise type error   
%
% Based on semantics specified for get_constraints_number for generic
% interface, but with type error if variable is not integer or variable

get_domain_size(I, Size) :-
        integer(I), !,
        Size = 1.
get_domain_size(V, Size) :-
        free(V), !,
        Size = 1.0Inf.
get_domain_size(_{gfd:Attr}, Size) ?- !,
        ( nonvar(Attr) ->
            Attr = gfd{prob:H, idx:Idx}, 
            restore_space_if_needed(H, SpH),
            g_get_var_domain_size(SpH, Idx, Size)
        ;
            Size = 1.0Inf
        ).
get_domain_size(X, Size) :-
        error(5, get_domain_size(X, Size)).


get_delta(I, Width) :-
        integer(I), !,
        Width = 0.
get_delta(V, Width) :-
        free(V), !,
        Width = 1.0Inf.
get_delta(_{gfd:Attr}, Width) ?- !,
        ( nonvar(Attr) ->
            Attr = gfd{prob:H, idx:Idx}, 
            restore_space_if_needed(H, SpH),
            g_get_var_domain_width(SpH, Idx, Width)
        ;
            Width = 1.0Inf
        ).
get_delta(V, W) :-
        error(5, get_delta(V, W)).


get_constraints_number(T, Size) :-
        integer(T), !,
        Size = 1.0Inf.
get_constraints_number(V, Size) :-
        free(V), !,
        Size = 0.
get_constraints_number(_{gfd:Attr}, Size) ?- !,
        ( nonvar(Attr) ->
            Attr = gfd{prob:H, idx:Idx},
            restore_space_if_needed(H, SpH),
            g_get_var_degree(SpH, Idx, Size)
        ;
            Size = 0
        ).
get_constraints_number(V, Size) :-
        error(5, get_constraints_number(V, Size)).


get_weighted_degree(T, Count) :-
        integer(T), !,
        Count = 1.0Inf.
get_weighted_degree(T, Count) :-
        free(T), !,
        Count = 0.
get_weighted_degree(_{gfd:Attr}, Count) ?- !,
        ( nonvar(Attr) ->
            Attr = gfd{prob:H, idx:Idx}, 
            restore_space_if_needed(H, SpH),
            g_get_var_afc(SpH, Idx, Count)
        ;
            Count = 0
        ).
get_weighted_degree(T, Count) :-
        error(5, get_weighted_degree(T, Count)).


get_regret_lwb(V, Count) :-
        integer(V), !,
        Count = 1.0Inf.
get_regret_lwb(V, Count) :-
        free(V), !,
        Count = 1.
get_regret_lwb(_{gfd:Attr}, Count) ?- !,
        ( nonvar(Attr) ->
            Attr = gfd{prob:H, idx:Idx},
            restore_space_if_needed(H, SpH),
            g_get_var_regret_lwb(SpH, Idx, Count)
        ;
            Count = 1
        ).
get_regret_lwb(V, Count) :-
        error(5, get_regret_lwb(V, Count)).


get_regret_upb(V, Count) :-
        nonvar(V), !,
        Count = 1.0Inf.
get_regret_upb(V, Count) :-
        free(V), !,
        Count = 1.
get_regret_upb(_{gfd:Attr}, Count) ?- !,
        ( nonvar(Attr) ->
            Attr = gfd{prob:H, idx:Idx},
            restore_space_if_needed(H, SpH),
            g_get_var_regret_upb(SpH, Idx, Count)
        ;
            Count = 1
        ).
get_regret_upb(V, Count) :-
        error(5, get_regret_upb(V, Count)).


impose_min(I, Min) :-
        integer(I), !,
        Min =< I.
impose_min(_{gfd:Attr}, Min0) ?- 
        nonvar(Attr),
        ( integer(Min0) -> Min = Min0
        ; number(Min0) ->  Min is integer(ceiling(Min0))
        ; fail
        ), !,
        Attr = gfd{prob:H, idx:Idx,bool:BI},
        gfdvar(Idx,BI, GV),
        post_new_event_no_wake(post_lwb(GV, Min), H).
impose_min(V, Min0) :-
        var(V),  
        ( integer(Min0) -> Min = Min0
        ; number(Min0) ->  Min is integer(ceiling(Min0))
        ; fail
        ), !,
        get_prob_handle_nvars(H, NV0),
        new_gfdvar(V, H, NV0,NV, _),
        setarg(nvars of gfd_prob, H, NV),
        gfd_default_interval(_, Max),
        post_new_event_no_wake(newvars_interval(NV,Min,Max), H).
impose_min(V, Min) :-
        error(5, impose_min(V, Min)).


impose_max(I, Max) :-
        integer(I), !,
        integer(Max),
        Max >= I.
impose_max(_{gfd:Attr}, Max0) ?- 
        nonvar(Attr), 
        ( integer(Max0) -> Max = Max0
        ; number(Max0) ->  Max is integer(floor(Max0))
        ; fail
        ), !,
        Attr = gfd{prob:H, idx:Idx,bool:BI},
        gfdvar(Idx,BI, GV),
        post_new_event_no_wake(post_upb(GV, Max), H).
impose_max(V, Max0) :-
        var(V),  !,
        ( integer(Max0) -> Max = Max0
        ; number(Max0) ->  Max is integer(floor(Max0))
        ; fail
        ), !,
        get_prob_handle_nvars(H, NV0),
        new_gfdvar(V, H, NV0,NV, _),
        setarg(nvars of gfd_prob, H, NV),
        gfd_default_interval(Min, _),
        post_new_event_no_wake(newvars_interval(NV,Min,Max), H).
impose_max(V, Max) :-
        error(5, impose_max(V, Max)).


impose_bounds(I, Min, Max) :-
        integer(I), !,
        Max >= I, Min =< I.
impose_bounds(V, Min0, Max0) :-
        var(V),
        ( integer(Min0) -> Min = Min0
        ; number(Min0) ->  Min is integer(ceiling(Min0))
        ; fail
        ), 
        ( integer(Max0) -> Max = Max0
        ; number(Max0) ->  Max is integer(floor(Max0))
        ; fail
        ), !,
        ( is_solver_var(V) ->
            gfd_set_var_bounds(V, Min, Max) 
            % wake % need explicit wake here (no longer!)
        ;
            get_prob_handle_nvars(H, N0),
            new_gfdvar(V, H, N0,N, _GV),
            do_update_newvars_with_domain_interval(H, N, Min, Max)
        ).
impose_bounds(V, Min, Max) :-
        error(5, impose_bounds(V, Min, Max)).


impose_domain(X, _Y{gfd:Attr}) ?- !,
        var(X),
        impose_domain1(X, Attr).
impose_domain(_X, Y) :-
        var(Y).
impose_domain(X, Y) :-
        nonvar(Y),
        X = Y.

    impose_domain1(_X, AttrY) :-
    	var(AttrY).
    impose_domain1(X, AttrY) :-
    	nonvar(AttrY),
        AttrY = gfd{idx:IdxY},
        get_prob_handle(H), 
        restore_space_if_needed(H, SpH),
        g_get_var_domain_handle(SpH, IdxY, Dom),
        impose_domain_on_var(X, H, Dom).

    impose_domain_on_var(_X{gfd:Attr}, H, Dom) ?-
        nonvar(Attr), !,
        Attr = gfd{idx:Idx, bool:BI},
        gfdvar(Idx,BI, GX),
        post_new_event(post_dom_handle([](GX), Dom), H).
    impose_domain_on_var(X, H, Dom) :-
        H = gfd_prob{nvars:N0},
        new_gfdvar(X, H, N0,N, _GX),
        setarg(nvars of gfd_prob, H, N),
        post_new_event(newvars_dom_handle(N,Dom), H).


exclude(I, Excl) ?-
        integer(I), !,
        (integer(Excl) -> true ; error(5, exclude(I, Excl))),
        I \== Excl.
exclude(V{gfd:Attr}, I) ?- 
        nonvar(Attr),
        (integer(I) -> true ; error(5, exclude(V, I))),
        Attr = gfd{prob:H, idx:Idx},
        post_new_event_no_wake(post_exclude_var_val(Idx, I), H).


exclude_range(V{gfd:Attr}, Lo, Hi) ?-
        nonvar(Attr), !,
        (integer(Hi) -> true ; error(5, exclude_range(V, Lo, Hi))),
        (integer(Lo) -> true ; error(5, exclude_range(V, Lo, Hi))),
        Attr = gfd{prob:H, idx:Idx,bool:BI},
        gfdvar(Idx,BI, GV),
        post_new_event_no_wake(post_var_interval_reif(GV, Lo, Hi, 0), H).
exclude_range(X, Lo, Hi) :- 
	integer(X),
	integer(Lo),
	integer(Hi),
	!,
	\+ ( Lo =< X, X =< Hi).
exclude_range(X, Lo, Hi) :-
	error(6, exclude_range(X, Lo, Hi)).



gfd_vars_impose_min(Vs, Min) :-
        gfd_default_interval(_, Max),
        impose_vars_bound(Vs, Min, Max, (#>=), Min).

gfd_vars_impose_max(Vs, Max) :-
        gfd_default_interval(Min, _),
        impose_vars_bound(Vs, Min, Max, (#=<), Max).

impose_vars_bound(Vs, Min, Max, Rel, Bound) :-
        collection_to_list(Vs, VList),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_oldvarlist_bounds1(VList, H, Min, Max, NV0, NV, GOList),
        (NV > NV0 ->
            setarg(nvars of gfd_prob, H, NV),
            post_new_event_no_wake(newvars_interval(NV,Min,Max), H)
        ;
            true
        ),
        ( GOList == [] ->
            true
        ;
            GOArr =.. [[]|GOList], 
            post_new_event_no_wake(post_collection_rel(default, GOArr, Rel, Bound), H)
        ).

gfd_vars_impose_bounds(Vs, Min, Max) :-
        collection_to_list(Vs, VList),
        get_prob_handle_nvars(H, NV0),
        ec_to_gecode_oldvarlist_bounds1(VList, H, Min, Max, NV0, NV, GOList),
        ( NV > NV0 ->
            % have new variables
            setarg(nvars of gfd_prob, H, NV),
            post_new_event_no_wake(newvars_interval(NV,Min,Max), H)
        ;
            true
        ),
        ( GOList == [] ->
            true
        ;
            GOArr =.. [[]|GOList], 
            post_new_event_no_wake(post_interval(GOArr, Min, Max), H)

        ).


:- tool(gfd_vars_impose_domain/2, gfd_vars_impose_domain/3).

gfd_vars_impose_domain(Vs, [Val], Module) ?-
        integer(Val), !,
        gfd_vars_impose_domain(Vs, Val, Module).
gfd_vars_impose_domain(Vs, Val, _Moudle) :-
        integer(Val), !,
        collection_to_list(Vs, VList),
        ( foreach(Val, VList), param(Val) do true ).
gfd_vars_impose_domain(Vs, Dom, Module) :-
        collection_to_list(Vs, VList),
        get_prob_handle_nvars(H, NV0),
        ( foreach(V, VList),
          fromto(NV0, NV1,NV2, NV),
          param(H),
          fromto(GVs, GVs1,GVs2, [])
        do
            ( integer(V) ->
                NV1 = NV2,
                GVs1 = [V|GVs2]
            ; var(V) ->
                ( get_gecode_attr(V, gfd{idx:I,bool:B}) ->
                    gfdvar(I, B, GV),
                    GVs1 = [GV|GVs2],
                    NV1 = NV2
                ; /* new gfd var */
                    new_gfdvar(V, H, NV1,NV2, _GV),
                    GVs1 = GVs2
                )
            ;
                fail
            )
        ),
        GVArr =.. [[]|GVs],
        impose_domain_on_vars(GVArr, Dom, H, NV0, NV, Module).



impose_domain_on_vars(GVArr, _{gfd:gfd{idx:Idx}}, H, NV0, NV, _Module) ?- !,
        restore_space_if_needed(H, SpH),
        g_get_var_domain_handle(SpH, Idx, Dom),
        (NV > NV0 ->
            setarg(nvars of gfd_prob, H, NV),
            post_new_event_no_wake(newvars_dom_handle(NV,Dom), H)
        ;
            true
        ),
        ( GVArr == [] ->
            true
        ;
            post_new_event_no_wake(post_dom_handle(GVArr, Dom), H)
        ).
impose_domain_on_vars(GVArr, Dom, H, NV0, NV, Module) :-
        compound(Dom), !,
        ( foreach(D,Dom), param(Module),
          foreach(Lo..Hi, Vals)
        do
            subdomain(D, Lo, Hi, Module)
        ),
        ValsArr =.. [[]|Vals],
        (NV > NV0 ->
            setarg(nvars of gfd_prob, H, NV),
            post_new_event_no_wake(newvars_dom(NV,ValsArr), H)
        ;
            true
        ),
        ( GVArr == [] ->
            true
        ;
            post_new_event_no_wake(post_dom(GVArr, ValsArr), H)
        ).


gfd_vars_exclude(Vars, Excl) :-
        (collection_to_list(Vars, VL) -> true ; error(5, gfd_vars_exclude(Vars, Excl))),
        (integer(Excl) -> true ; error(5, gfd_vars_exclude(Vars, Excl))),
        ( foreach(V, VL),
          param(Vars, Excl),
          fromto(GVars, GVs1,GVs2, [])
        do
            ( integer(V) -> 
                V \== Excl,
                GVs1 = GVs2
            ; get_gecode_attr(V, Attr) -> 
                Attr = gfd{idx:Idx},
                GVs1 = [Idx|GVs2]
            ;
                 error(6, gfd_vars_exclude(Vars, Excl))
            )
        ),
        GVArr =.. [[]|GVars],
        get_prob_handle(H),
        post_new_event_no_wake(post_exclude_val(GVArr, Excl), H).


gfd_vars_exclude_range(Vars, Lo, Hi) :-
        collection_to_list(Vars, VL),
        integer(Hi),
        integer(Lo),
        !,
        Hi >= Lo,
        ( foreach(V, VL),
          param(Vars, Lo, Hi),
          fromto(GVars, GVs1,GVs2, [])
        do
            ( integer(V) -> 
                \+ ( Lo =< V, V =< Hi),
                GVs1 = GVs2
            ; get_gecode_attr(V, Attr) -> 
                Attr = gfd{idx:I},
                GVs1 = [I|GVs2]
            ;
                 error(6, gfd_vars_exclude_range(Vars, Lo, Hi))
            )
        ),
        get_prob_handle(H),
        GVArr =.. [[]|GVars],
        post_new_event_no_wake(post_exclude_range(GVArr, Lo, Hi), H).
gfd_vars_exclude_range(Vars, Lo, Hi) :-
        error(6, gfd_vars_exclude_range(Vars, Lo, Hi)).


:- tool(gfd_vars_exclude_domain/2, gfd_vars_exclude_domain/3).

gfd_vars_exclude_domain(Vars, Dom, Module) :-
        check_collection_to_list(Vars, VL), 
        ( foreach(V, VL),
          foreach(GV, GVars)
        do
            get_gecode_var_or_int(V, GV) 
        ), !,
        (GVars == [] -> 
            true
        ;
            GVArr =.. [[]|GVars],
            exclude_domain1(GVArr, Dom, Module)
        ).
gfd_vars_exclude_domain(Vs, Dom, _M) :-
        ( get_bip_error(E) -> true ; E = 6),
        error(E, exclude_domain(Vs, Dom)).

exclude_domain1(Vs, _{gfd:gfd{idx:Idx}}, _Module) ?- !,
        get_prob_handle(H),
        restore_space_if_needed(H, SpH),
        g_get_var_domain_handle(SpH, Idx, Dom),
        post_new_event_no_wake(post_exclude_dom_handle(Vs, Dom), H).
exclude_domain1(Vs, Dom, Module) :-
        compound(Dom),
        ( foreach(D,Dom), param(Module),
          foreach(Lo..Hi, Vals)
        do
            subdomain(D, Lo, Hi, Module)
        ),
        ValsArr =.. [[]|Vals],
        get_prob_handle(H),
        post_new_event_no_wake(post_exclude_dom(Vs, ValsArr), H).


is_solver_type(I) :- integer(I), !.
is_solver_type(_{gfd:Attr}) ?- 
        nonvar(Attr),
        Attr = gfd{}.

is_solver_var(_{gfd:Attr}) ?-
        nonvar(Attr),
        Attr = gfd{}.

is_exact_solver_var(V) :-
        is_solver_var(V).

integers(V) :-
        get_prob_handle(H),
        ( var(V) ->
            ec_to_gecode_var(V, H, _)
        ;
            collection_to_list(V, VList),
            ec_to_gecode_varlist(VList, H, _, _)
        ).


msg(X, Y, Dom) :-
        ( integer(X), integer(Y), X == Y ->
            % X =:= Y, Msg is a singleton, so directly bind it
            Dom = X
        ;
            ( ec_to_gecode_oldvar(X, GX),
              ec_to_gecode_oldvar(Y, GY) ->
                get_prob_handle(H),
                restore_space_if_needed(H, _),
                H = gfd_prob{nvars:N0},
                % The new domain must be added to a new gfd var
                % use Dom0 in case Dom is an exising Domain var
                new_gfdvar(Dom0, H, N0,N, _GDom0),
                setarg(nvars of gfd_prob, H, N),
                % we follow Gecode's recomputation style, the
                % union will be recomputed, rather than storing it 
                post_new_event(newvars_dom_union(GX,GY,N), H),
                Dom0 = Dom
            ;
                true
            )
        ).


is_in_domain(Val, Var) :-
        integer(Val), !,
        ( var(Var) ->
            get_gecode_attr(Var, Attr),
            Attr = gfd{idx:Idx,prob:H},
            restore_space_if_needed(H, SpH),
            g_check_val_is_in_var_domain(SpH, Idx, Val)
        ;
            Val == Var
        ).

is_in_domain(Val, Var, Result) :-
        (is_in_domain(Val, Var) -> Result = yes ; Result = no).


% note these don't need a space, so no need to update handle
gfd_maxint(X) :-
        g_get_gfd_maxint(X).

gfd_minint(X) :-
        g_get_gfd_minint(X).


:- erase_module(gfd_gac),
   create_constraint_pool(gfd_gac, 0, [
      (#\=)/2 -> '#\\=_c'/3,
      (#=)/2 -> '#=_c'/3,
      (#<)/2 -> '#<_c'/3,
      (#>)/2 -> '#>_c'/3,
      (#>=)/2 -> '#>=_c'/3,
      (#=<)/2 -> '#=<_c'/3,
      (and)/2 -> and_c/3,
      (or)/2 -> or_c/3,
      (xor)/2 -> xor_c/3,
      neg/1 -> neg_c/2,
      '<=>'/2 -> '<=>_c'/3,
      '=>'/2 -> '=>_c'/3,
/* currently gecode does not support gac for reified expressions    
      (#\=)/3 -> '#\\=_reif_c'/4,
      (#=)/3 -> '#=_reif_c'/4,
      (#<)/3 -> '#<_reif_c'/4,
      (#>)/3 -> '#>_reif_c'/4,
      (#>=)/3 -> '#>=_reif_c'/4,
      (#=<)/3 -> '#=<_reif_c'/4,
      (and)/3 -> and_reif_c/4,
      (or)/3 -> or_reif_c/4,
      (xor)/3 -> xor_reif_c/4,
      neg/2 -> neg_c/3,
      '<=>'/3 -> '<=>_reif_c'/4,
      '=>'/3 -> '=>_reif_c'/4, */
      ordered/2 -> ordered_c/3,
      among/4 -> among_c/5,
      count_matches/4 -> count_matches_c/5,
      mem/2 -> mem_c/3,
      mem/3 -> mem_reif_c/4,
      all_le/2 -> all_le_c/3,
      all_lt/2 -> all_lt_c/3,
      all_ge/2 -> all_ge_c/3,
      all_gt/2 -> all_gt_c/3,
      all_eq/2 -> all_eq_c/3,
      all_ne/2 -> all_ne_c/3,
      lex_le/2 -> lex_le_c/3,
      lex_lt/2 -> lex_lt_c/3,
      lex_ge/2 -> lex_ge_c/3,
      lex_gt/2 -> lex_gt_c/3,
      lex_eq/2 -> lex_eq_c/3,
      lex_ne/2 -> lex_ne_c/3,
      alldifferent/1 -> alldifferent_c/2,
      alldifferent_cst/2 -> alldifferent_cst_c/3,
      bool_channeling/3 -> bool_channeling_c/4,
      count/4 -> count_c/5,
      occurrences/3 -> occurrences_c/4,
      atmost/3 -> atmost_c/4,
      atleast/3 -> atleast_c/4,
      element/3 -> element_c/4,
      element_g/3 -> element_g_c/4,
      circuit/1 -> circuit_c/2,
      circuit/3 -> circuit_c/4,
      circuit/4 -> circuit_c/5,
      circuit_g/1 -> circuit_g_c/2,
      circuit_g/3 -> circuit_g_c/4,
      circuit_g/4 -> circuit_g_c/5,
      circuit_offset/2 -> circuit_offset_c/3,
      circuit_offset/4 -> circuit_offset_c/5,
      circuit_offset/5 -> circuit_offset_c/6,
      circuit_offset_g/2 -> circuit_offset_g_c/3,
      circuit_offset_g/4 -> circuit_offset_g_c/5,
      circuit_offset_g/5 -> circuit_offset_g_c/6,
      ham_path/3 -> ham_path_c/4,
      ham_path/5 -> ham_path_c/6,
      ham_path/6 -> ham_path_c/7,
      ham_path_g/3 -> ham_path_g_c/4,
      ham_path_g/5 -> ham_path_g_c/6,
      ham_path_g/6 -> ham_path_g_c/7,
      ham_path_offset/4 -> ham_path_offset_c/5,
      ham_path_offset/6 -> ham_path_offset_c/7,
      ham_path_offset/7 -> ham_path_offset_c/8,
      ham_path_offset_g/4 -> ham_path_offset_g_c/5,
      ham_path_offset_g/6 -> ham_path_offset_g_c/7,
      ham_path_offset_g/7 -> ham_path_offset_g_c/8,
      gcc/2 -> gcc_c/3,
      precede/3 -> precede_c/4,
      inverse/2 -> inverse_c/3,
      inverse_g/2 -> inverse_g_c/3,
      inverse/4 -> inverse_c/5,
      inverse_g/4 -> inverse_g_c/5,
      min/2 -> minlist_c/3,
      max/2 -> maxlist_c/3,
      sumlist/2 -> sum_c/3,
      sum/2 -> sum_c/3,
      sum/3 -> sum_c/4,
      scalar_product/4 -> scalar_product_c/5,
      sequence/4 -> sequence_c/5,
      sequence/5 -> sequence_c/6,
      extensional/4 -> extensional_c/5,
      regular/2 ->regular_c/3,
      table/2 -> table_c/3,
      table/3 -> table_c/4
   ]).

:- erase_module(gfd_bc),
   create_constraint_pool(gfd_bc, 0, [
      (#\=)/2 -> '#\\=_c'/3,
      (#=)/2 -> '#=_c'/3,
      (#<)/2 -> '#<_c'/3,
      (#>)/2 -> '#>_c'/3,
      (#>=)/2 -> '#>=_c'/3,
      (#=<)/2 -> '#=<_c'/3,
      (and)/2 -> and_c/3,
      (or)/2 -> or_c/3,
      (xor)/2 -> xor_c/3,
      neg/1 -> neg_c/2,
      '<=>'/2 -> '<=>_c'/3,
      '=>'/2 -> '=>_c'/3,
      (#\=)/3 -> '#\\=_reif_c'/4,
      (#=)/3 -> '#=_reif_c'/4,
      (#<)/3 -> '#<_reif_c'/4,
      (#>)/3 -> '#>_reif_c'/4,
      (#>=)/3 -> '#>=_reif_c'/4,
      (#=<)/3 -> '#=<_reif_c'/4,
      (and)/3 -> and_reif_c/4,
      (or)/3 -> or_reif_c/4,
      (xor)/3 -> xor_reif_c/4,
      neg/2 -> neg_c/3,
      '<=>'/3 -> '<=>_reif_c'/4,
      '=>'/3 -> '=>_reif_c'/4,
      divmod/4 -> divmod_c/5,
      all_le/2 -> all_le_c/3,
      all_lt/2 -> all_lt_c/3,
      all_ge/2 -> all_ge_c/3,
      all_gt/2 -> all_gt_c/3,
      all_eq/2 -> all_eq_c/3,
      all_ne/2 -> all_ne_c/3,
      ordered/2 -> ordered_c/3,
      lex_le/2 -> lex_le_c/3,
      lex_lt/2 -> lex_lt_c/3,
      lex_ge/2 -> lex_ge_c/3,
      lex_gt/2 -> lex_gt_c/3,
      lex_eq/2 -> lex_eq_c/3,
      lex_ne/2 -> lex_ne_c/3,
      alldifferent/1 -> alldifferent_c/2,
      alldifferent_cst/2 -> alldifferent_cst_c/3,
      element/3 -> element_c/4,
      element_g/3 -> element_g_c/4,
      gcc/2 -> gcc_c/3,
      count/4 -> count_c/5,
      occurrences/3 -> occurrences_c/4,
      min/2 -> minlist_c/3,
      max/2 -> maxlist_c/3,
      sumlist/2 -> sum_c/3,
      sum/2 -> sum_c/3,
      sum/3 -> sum_c/4,
      sum/4 -> sum_reif_c/5,
      scalar_product/4 -> scalar_product_c/5,
      scalar_product/5 -> scalar_product_reif_c/6,
      sorted_g/3 -> sorted_g_c/4,
      sorted/2 -> sorted_c/3,
      sorted/3 -> sorted_c/4
  ]).

:- erase_module(gfd_vc),
   create_constraint_pool(gfd_vc, 0, [
      alldifferent/1 -> alldifferent_c/2,
      alldifferent_cst/2 -> alldifferent_cst_c/3,
      circuit/1 -> circuit_c/2,
      circuit/3 -> circuit_c/4,
      circuit/4 -> circuit_c/5,
      circuit_g/1 -> circuit_g_c/2,
      circuit_g/3 -> circuit_g_c/4,
      circuit_g/4 -> circuit_g_c/5,
      circuit_offset/2 -> circuit_offset_c/3,
      circuit_offset/4 -> circuit_offset_c/5,
      circuit_offset/5 -> circuit_offset_c/6,
      circuit_offset_g/2 -> circuit_offset_g_c/3,
      circuit_offset_g/4 -> circuit_offset_g_c/5,
      circuit_offset_g/5 -> circuit_offset_g_c/6,
      ham_path/3 -> ham_path_c/4,
      ham_path/5 -> ham_path_c/6,
      ham_path/6 -> ham_path_c/7,
      ham_path_g/3 -> ham_path_g_c/4,
      ham_path_g/5 -> ham_path_g_c/6,
      ham_path_g/6 -> ham_path_g_c/7,
      ham_path_offset/4 -> ham_path_offset_c/5,
      ham_path_offset/6 -> ham_path_offset_c/7,
      ham_path_offset/7 -> ham_path_offset_c/8,
      ham_path_offset_g/4 -> ham_path_offset_g_c/5,
      ham_path_offset_g/6 -> ham_path_offset_g_c/7,
      ham_path_offset_g/7 -> ham_path_offset_g_c/8,
%      cumulative/4 -> cumulative_c/5,
      cumulatives/5 -> cumulatives_c/6,
      cumulatives_min/5 -> cumulatives_min_c/6,
      cumulatives_g/5 -> cumulatives_g_c/6,
      cumulatives_min_g/5 -> cumulatives_min_g_c/6,
      inverse/2 -> inverse_c/3,
      inverse_g/2 -> inverse_g_c/3,
      inverse/4 -> inverse_c/5,
      inverse_g/4 -> inverse_g_c/5,
      gcc/2 -> gcc_c/3,
      ordered/2 -> ordered_c/3
  ]).


:- comment(include, "gfd_comments.ecl").
