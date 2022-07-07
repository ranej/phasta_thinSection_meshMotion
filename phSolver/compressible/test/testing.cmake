set(casename staticBubble)
set(CDIR ${CASES}/${casename}/run)
set(NPROCS 8)
set(tol 1e-6)
c_serial_test(inpCfg_${casename} ${CDIR} cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})
c_serial_test(resCfg_${casename} ${CDIR} cp -r ${CDIR}/../chef/${NPROCS}-procs_case ${CDIR})
c_parallel_test(${casename} ${NPROCS} ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
set(compareArgs
  ${CDIR}/${NPROCS}-procs_case/
  ${CDIR}/${NPROCS}-procs_case_ref/
  ${NPROCS} ${tol})
c_parallel_test(compareRestart-${casename} ${NPROCS} ${CDIR}
	${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})
