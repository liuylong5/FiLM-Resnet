! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The ODE Jacobian of Chemical Model File
! 
! Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : aqchem_Jacobian.F90
! Time                 : Tue Dec  2 14:12:20 2014
! Working directory    : /home/kfahey/kpp-2.2.3/ASSUMPTION_TEST/CLEAN_INIT/TESTS/CMAQ_BASE_022414/3D/FINAL_120114
! Equation file        : aqchem.kpp
! Output root filename : aqchem
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE aqchem_Jacobian

  USE aqchem_Parameters
  USE aqchem_JacobianSP

  IMPLICIT NONE

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Jac_SP - the Jacobian of Variables in sparse matrix representation
!   Arguments :
!      V         - Concentrations of variable species (local)
!      F         - Concentrations of fixed species (local)
!      RCT       - Rate constants (local)
!      JVS       - sparse Jacobian of variables
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Jac_SP ( V, F, RCT, JVS )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! F - Concentrations of fixed species (local)
  REAL(kind=dp) :: F(NFIX)
! RCT - Rate constants (local)
  REAL(kind=dp) :: RCT(NREACT)
! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)


! Local variables
! B - Temporary array
  REAL(kind=dp) :: B(142)

! B(1) = dA(1)/dV(59)
  B(1) = RCT(1)
! B(2) = dA(2)/dV(46)
  B(2) = RCT(2)
! B(3) = dA(3)/dV(47)
  B(3) = RCT(3)
! B(4) = dA(4)/dV(48)
  B(4) = RCT(4)
! B(5) = dA(5)/dV(49)
  B(5) = RCT(5)
! B(6) = dA(6)/dV(50)
  B(6) = RCT(6)
! B(7) = dA(7)/dV(51)
  B(7) = RCT(7)
! B(8) = dA(8)/dV(52)
  B(8) = RCT(8)
! B(9) = dA(9)/dV(57)
  B(9) = RCT(9)
! B(10) = dA(10)/dV(58)
  B(10) = RCT(10)
! B(11) = dA(11)/dV(54)
  B(11) = RCT(11)
! B(12) = dA(12)/dV(56)
  B(12) = RCT(12)
! B(13) = dA(13)/dV(82)
  B(13) = RCT(13)
! B(14) = dA(14)/dV(60)
  B(14) = RCT(14)
! B(15) = dA(15)/dV(71)
  B(15) = RCT(15)
! B(16) = dA(16)/dV(62)
  B(16) = RCT(16)
! B(17) = dA(17)/dV(74)
  B(17) = RCT(17)
! B(18) = dA(18)/dV(80)
  B(18) = RCT(18)
! B(19) = dA(19)/dV(65)
  B(19) = RCT(19)
! B(20) = dA(20)/dV(77)
  B(20) = RCT(20)
! B(21) = dA(21)/dV(75)
  B(21) = RCT(21)
! B(22) = dA(22)/dV(68)
  B(22) = RCT(22)
! B(23) = dA(23)/dV(53)
  B(23) = RCT(23)
! B(24) = dA(24)/dV(55)
  B(24) = RCT(24)
! B(25) = dA(25)/dV(1)
  B(25) = RCT(25)
! B(26) = dA(26)/dV(2)
  B(26) = RCT(26)
! B(27) = dA(27)/dV(3)
  B(27) = RCT(27)
! B(28) = dA(28)/dV(4)
  B(28) = RCT(28)
! B(29) = dA(29)/dV(5)
  B(29) = RCT(29)
! B(30) = dA(30)/dV(6)
  B(30) = RCT(30)
! B(31) = dA(31)/dV(7)
  B(31) = RCT(31)
! B(32) = dA(32)/dV(8)
  B(32) = RCT(32)
! B(33) = dA(33)/dV(82)
  B(33) = RCT(33)
! B(34) = dA(34)/dV(81)
  B(34) = RCT(34)
! B(35) = dA(35)/dV(60)
  B(35) = RCT(35)
! B(36) = dA(36)/dV(71)
  B(36) = RCT(36)
! B(37) = dA(37)/dV(72)
  B(37) = RCT(37)
! B(38) = dA(38)/dV(62)
  B(38) = RCT(38)
! B(39) = dA(39)/dV(65)
  B(39) = RCT(39)
! B(40) = dA(40)/dV(68)
  B(40) = RCT(40)
! B(42) = dA(42)/dV(67)
  B(42) = RCT(42)
! B(43) = dA(43)/dV(73)
  B(43) = RCT(43)
! B(44) = dA(44)/dV(79)
  B(44) = RCT(44)*V(81)
! B(45) = dA(44)/dV(81)
  B(45) = RCT(44)*V(79)
! B(46) = dA(45)/dV(78)
  B(46) = RCT(45)*V(79)
! B(47) = dA(45)/dV(79)
  B(47) = RCT(45)*V(78)
! B(48) = dA(46)/dV(61)
  B(48) = RCT(46)*V(79)
! B(49) = dA(46)/dV(79)
  B(49) = RCT(46)*V(61)
! B(50) = dA(47)/dV(72)
  B(50) = RCT(47)*V(79)
! B(51) = dA(47)/dV(79)
  B(51) = RCT(47)*V(72)
! B(52) = dA(48)/dV(70)
  B(52) = RCT(48)*V(79)
! B(53) = dA(48)/dV(79)
  B(53) = RCT(48)*V(70)
! B(54) = dA(49)/dV(63)
  B(54) = RCT(49)*V(64)
! B(55) = dA(49)/dV(64)
  B(55) = RCT(49)*V(63)
! B(56) = dA(50)/dV(66)
  B(56) = RCT(50)*V(79)
! B(57) = dA(50)/dV(79)
  B(57) = RCT(50)*V(66)
! B(58) = dA(51)/dV(69)
  B(58) = RCT(51)*V(79)
! B(59) = dA(51)/dV(79)
  B(59) = RCT(51)*V(69)
! B(60) = dA(52)/dV(64)
  B(60) = RCT(52)*V(79)
! B(61) = dA(52)/dV(79)
  B(61) = RCT(52)*V(64)
! B(62) = dA(53)/dV(73)
  B(62) = RCT(53)*V(79)
! B(63) = dA(53)/dV(79)
  B(63) = RCT(53)*V(73)
! B(64) = dA(54)/dV(76)
  B(64) = RCT(54)*V(79)
! B(65) = dA(54)/dV(79)
  B(65) = RCT(54)*V(76)
! B(66) = dA(55)/dV(53)
  B(66) = RCT(55)*F(2)
! B(68) = dA(56)/dV(55)
  B(68) = RCT(56)*F(2)
! B(70) = dA(57)/dV(74)
  B(70) = RCT(57)*V(79)*V(81)
! B(71) = dA(57)/dV(79)
  B(71) = RCT(57)*V(74)*V(81)
! B(72) = dA(57)/dV(81)
  B(72) = RCT(57)*V(74)*V(79)
! B(73) = dA(58)/dV(80)
  B(73) = RCT(58)*V(82)
! B(74) = dA(58)/dV(82)
  B(74) = RCT(58)*V(80)
! B(75) = dA(59)/dV(80)
  B(75) = RCT(59)*V(81)
! B(76) = dA(59)/dV(81)
  B(76) = RCT(59)*V(80)
! B(77) = dA(60)/dV(78)
  B(77) = RCT(60)*V(80)
! B(78) = dA(60)/dV(80)
  B(78) = RCT(60)*V(78)
! B(79) = dA(61)/dV(13)
  B(79) = RCT(61)*V(82)
! B(80) = dA(61)/dV(82)
  B(80) = RCT(61)*V(13)
! B(81) = dA(62)/dV(13)
  B(81) = RCT(62)*V(81)
! B(82) = dA(62)/dV(81)
  B(82) = RCT(62)*V(13)
! B(83) = dA(63)/dV(13)
  B(83) = RCT(63)*V(78)
! B(84) = dA(63)/dV(78)
  B(84) = RCT(63)*V(13)
! B(85) = dA(64)/dV(12)
  B(85) = RCT(64)*V(82)
! B(86) = dA(64)/dV(82)
  B(86) = RCT(64)*V(12)
! B(87) = dA(65)/dV(12)
  B(87) = RCT(65)*V(81)
! B(88) = dA(65)/dV(81)
  B(88) = RCT(65)*V(12)
! B(89) = dA(66)/dV(12)
  B(89) = RCT(66)*V(78)
! B(90) = dA(66)/dV(78)
  B(90) = RCT(66)*V(12)
! B(91) = dA(67)/dV(12)
  B(91) = RCT(67)*V(13)*V(82)
! B(92) = dA(67)/dV(13)
  B(92) = RCT(67)*V(12)*V(82)
! B(93) = dA(67)/dV(82)
  B(93) = RCT(67)*V(12)*V(13)
! B(94) = dA(68)/dV(12)
  B(94) = RCT(68)*V(13)*V(81)
! B(95) = dA(68)/dV(13)
  B(95) = RCT(68)*V(12)*V(81)
! B(96) = dA(68)/dV(81)
  B(96) = RCT(68)*V(12)*V(13)
! B(97) = dA(69)/dV(12)
  B(97) = RCT(69)*V(13)*V(78)
! B(98) = dA(69)/dV(13)
  B(98) = RCT(69)*V(12)*V(78)
! B(99) = dA(69)/dV(78)
  B(99) = RCT(69)*V(12)*V(13)
! B(100) = dA(70)/dV(77)
  B(100) = RCT(70)*V(79)*V(81)
! B(101) = dA(70)/dV(79)
  B(101) = RCT(70)*V(77)*V(81)
! B(102) = dA(70)/dV(81)
  B(102) = RCT(70)*V(77)*V(79)
! B(103) = dA(71)/dV(75)
  B(103) = RCT(71)*V(79)*V(81)
! B(104) = dA(71)/dV(79)
  B(104) = RCT(71)*V(75)*V(81)
! B(105) = dA(71)/dV(81)
  B(105) = RCT(71)*V(75)*V(79)
! B(106) = dA(72)/dV(75)
  B(106) = RCT(72)*V(81)
! B(107) = dA(72)/dV(81)
  B(107) = RCT(72)*V(75)
! B(108) = dA(73)/dV(82)
  B(108) = RCT(73)
! B(109) = dA(74)/dV(60)
  B(109) = RCT(74)
! B(110) = dA(75)/dV(71)
  B(110) = RCT(75)
! B(111) = dA(76)/dV(62)
  B(111) = RCT(76)
! B(112) = dA(77)/dV(74)
  B(112) = RCT(77)
! B(113) = dA(78)/dV(80)
  B(113) = RCT(78)
! B(114) = dA(79)/dV(65)
  B(114) = RCT(79)
! B(115) = dA(80)/dV(77)
  B(115) = RCT(80)
! B(116) = dA(81)/dV(75)
  B(116) = RCT(81)
! B(117) = dA(82)/dV(67)
  B(117) = RCT(82)
! B(118) = dA(83)/dV(68)
  B(118) = RCT(83)
! B(119) = dA(84)/dV(53)
  B(119) = RCT(84)
! B(120) = dA(85)/dV(55)
  B(120) = RCT(85)
! B(121) = dA(86)/dV(76)
  B(121) = RCT(86)
! B(122) = dA(87)/dV(61)
  B(122) = RCT(87)
! B(123) = dA(88)/dV(63)
  B(123) = RCT(88)
! B(124) = dA(89)/dV(69)
  B(124) = RCT(89)
! B(125) = dA(90)/dV(9)
  B(125) = RCT(90)
! B(126) = dA(91)/dV(10)
  B(126) = RCT(91)
! B(127) = dA(92)/dV(39)
  B(127) = RCT(92)
! B(128) = dA(93)/dV(11)
  B(128) = RCT(93)
! B(129) = dA(94)/dV(14)
  B(129) = RCT(94)
! B(130) = dA(95)/dV(15)
  B(130) = RCT(95)
! B(131) = dA(96)/dV(16)
  B(131) = RCT(96)
! B(132) = dA(97)/dV(17)
  B(132) = RCT(97)
! B(133) = dA(98)/dV(12)
  B(133) = RCT(98)
! B(134) = dA(99)/dV(13)
  B(134) = RCT(99)
! B(135) = dA(100)/dV(79)
  B(135) = RCT(100)
! B(136) = dA(101)/dV(64)
  B(136) = RCT(101)
! B(137) = dA(102)/dV(81)
  B(137) = RCT(102)
! B(138) = dA(103)/dV(78)
  B(138) = RCT(103)
! B(139) = dA(104)/dV(72)
  B(139) = RCT(104)
! B(140) = dA(105)/dV(70)
  B(140) = RCT(105)
! B(141) = dA(106)/dV(66)
  B(141) = RCT(106)
! B(142) = dA(107)/dV(73)
  B(142) = RCT(107)

! Construct the Jacobian terms from B's
! JVS(1) = Jac_FULL(1,1)
  JVS(1) = -B(25)
! JVS(2) = Jac_FULL(2,2)
  JVS(2) = -B(26)
! JVS(3) = Jac_FULL(3,3)
  JVS(3) = -B(27)
! JVS(4) = Jac_FULL(4,4)
  JVS(4) = -B(28)
! JVS(5) = Jac_FULL(5,5)
  JVS(5) = -B(29)
! JVS(6) = Jac_FULL(6,6)
  JVS(6) = -B(30)
! JVS(7) = Jac_FULL(7,7)
  JVS(7) = -B(31)
! JVS(8) = Jac_FULL(8,8)
  JVS(8) = -B(32)
! JVS(9) = Jac_FULL(9,8)
  JVS(9) = B(32)
! JVS(10) = Jac_FULL(9,9)
  JVS(10) = -B(125)
! JVS(11) = Jac_FULL(10,6)
  JVS(11) = B(30)
! JVS(12) = Jac_FULL(10,10)
  JVS(12) = -B(126)
! JVS(13) = Jac_FULL(11,7)
  JVS(13) = B(31)
! JVS(14) = Jac_FULL(11,11)
  JVS(14) = -B(128)
! JVS(15) = Jac_FULL(12,12)
  JVS(15) = -B(133)
! JVS(16) = Jac_FULL(13,13)
  JVS(16) = -B(134)
! JVS(17) = Jac_FULL(14,14)
  JVS(17) = -B(129)
! JVS(18) = Jac_FULL(15,15)
  JVS(18) = -B(130)
! JVS(19) = Jac_FULL(16,16)
  JVS(19) = -B(131)
! JVS(20) = Jac_FULL(17,4)
  JVS(20) = B(28)
! JVS(21) = Jac_FULL(17,17)
  JVS(21) = -B(132)
! JVS(22) = Jac_FULL(18,18)
  JVS(22) = 0
! JVS(23) = Jac_FULL(18,78)
  JVS(23) = B(138)
! JVS(24) = Jac_FULL(18,81)
  JVS(24) = B(137)
! JVS(25) = Jac_FULL(18,82)
  JVS(25) = B(108)
! JVS(26) = Jac_FULL(19,19)
  JVS(26) = 0
! JVS(27) = Jac_FULL(19,60)
  JVS(27) = B(109)
! JVS(28) = Jac_FULL(20,20)
  JVS(28) = 0
! JVS(29) = Jac_FULL(20,70)
  JVS(29) = B(140)
! JVS(30) = Jac_FULL(20,71)
  JVS(30) = B(110)
! JVS(31) = Jac_FULL(20,72)
  JVS(31) = B(139)
! JVS(32) = Jac_FULL(21,21)
  JVS(32) = 0
! JVS(33) = Jac_FULL(21,62)
  JVS(33) = B(111)
! JVS(34) = Jac_FULL(22,22)
  JVS(34) = 0
! JVS(35) = Jac_FULL(22,74)
  JVS(35) = B(112)
! JVS(36) = Jac_FULL(23,23)
  JVS(36) = 0
! JVS(37) = Jac_FULL(23,80)
  JVS(37) = B(113)
! JVS(38) = Jac_FULL(24,24)
  JVS(38) = 0
! JVS(39) = Jac_FULL(24,65)
  JVS(39) = B(114)
! JVS(40) = Jac_FULL(24,66)
  JVS(40) = B(141)
! JVS(41) = Jac_FULL(25,25)
  JVS(41) = 0
! JVS(42) = Jac_FULL(25,77)
  JVS(42) = B(115)
! JVS(43) = Jac_FULL(26,26)
  JVS(43) = 0
! JVS(44) = Jac_FULL(26,75)
  JVS(44) = B(116)
! JVS(45) = Jac_FULL(27,27)
  JVS(45) = 0
! JVS(46) = Jac_FULL(27,67)
  JVS(46) = B(117)
! JVS(47) = Jac_FULL(27,73)
  JVS(47) = B(142)
! JVS(48) = Jac_FULL(27,76)
  JVS(48) = B(121)
! JVS(49) = Jac_FULL(28,28)
  JVS(49) = 0
! JVS(50) = Jac_FULL(28,68)
  JVS(50) = B(118)
! JVS(51) = Jac_FULL(29,29)
  JVS(51) = 0
! JVS(52) = Jac_FULL(29,53)
  JVS(52) = B(119)
! JVS(53) = Jac_FULL(30,30)
  JVS(53) = 0
! JVS(54) = Jac_FULL(30,55)
  JVS(54) = B(120)
! JVS(55) = Jac_FULL(31,31)
  JVS(55) = 0
! JVS(56) = Jac_FULL(31,61)
  JVS(56) = B(122)
! JVS(57) = Jac_FULL(32,32)
  JVS(57) = 0
! JVS(58) = Jac_FULL(32,63)
  JVS(58) = B(123)
! JVS(59) = Jac_FULL(33,33)
  JVS(59) = 0
! JVS(60) = Jac_FULL(33,69)
  JVS(60) = B(124)
! JVS(61) = Jac_FULL(34,9)
  JVS(61) = B(125)
! JVS(62) = Jac_FULL(34,34)
  JVS(62) = 0
! JVS(63) = Jac_FULL(35,12)
  JVS(63) = B(133)
! JVS(64) = Jac_FULL(35,35)
  JVS(64) = 0
! JVS(65) = Jac_FULL(36,13)
  JVS(65) = B(134)
! JVS(66) = Jac_FULL(36,36)
  JVS(66) = 0
! JVS(67) = Jac_FULL(37,10)
  JVS(67) = B(126)
! JVS(68) = Jac_FULL(37,37)
  JVS(68) = 0
! JVS(69) = Jac_FULL(38,38)
  JVS(69) = 0
! JVS(70) = Jac_FULL(38,39)
  JVS(70) = B(127)
! JVS(71) = Jac_FULL(39,39)
  JVS(71) = -B(127)
! JVS(72) = Jac_FULL(39,53)
  JVS(72) = 0.04*B(66)
! JVS(73) = Jac_FULL(39,55)
  JVS(73) = 0.04*B(68)
! JVS(74) = Jac_FULL(40,11)
  JVS(74) = B(128)
! JVS(75) = Jac_FULL(40,40)
  JVS(75) = 0
! JVS(76) = Jac_FULL(41,41)
  JVS(76) = 0
! JVS(77) = Jac_FULL(41,79)
  JVS(77) = B(135)
! JVS(78) = Jac_FULL(42,14)
  JVS(78) = B(129)
! JVS(79) = Jac_FULL(42,42)
  JVS(79) = 0
! JVS(80) = Jac_FULL(43,15)
  JVS(80) = B(130)
! JVS(81) = Jac_FULL(43,43)
  JVS(81) = 0
! JVS(82) = Jac_FULL(44,16)
  JVS(82) = B(131)
! JVS(83) = Jac_FULL(44,44)
  JVS(83) = 0
! JVS(84) = Jac_FULL(45,17)
  JVS(84) = B(132)
! JVS(85) = Jac_FULL(45,45)
  JVS(85) = 0
! JVS(86) = Jac_FULL(46,46)
  JVS(86) = -B(2)
! JVS(87) = Jac_FULL(46,60)
  JVS(87) = B(14)
! JVS(88) = Jac_FULL(47,47)
  JVS(88) = -B(3)
! JVS(89) = Jac_FULL(47,71)
  JVS(89) = B(15)
! JVS(90) = Jac_FULL(48,48)
  JVS(90) = -B(4)
! JVS(91) = Jac_FULL(48,62)
  JVS(91) = B(16)
! JVS(92) = Jac_FULL(49,49)
  JVS(92) = -B(5)
! JVS(93) = Jac_FULL(49,74)
  JVS(93) = B(17)
! JVS(94) = Jac_FULL(50,50)
  JVS(94) = -B(6)
! JVS(95) = Jac_FULL(50,80)
  JVS(95) = B(18)
! JVS(96) = Jac_FULL(51,51)
  JVS(96) = -B(7)
! JVS(97) = Jac_FULL(51,65)
  JVS(97) = B(19)
! JVS(98) = Jac_FULL(52,52)
  JVS(98) = -B(8)
! JVS(99) = Jac_FULL(52,77)
  JVS(99) = B(20)
! JVS(100) = Jac_FULL(53,53)
  JVS(100) = -B(23)-B(66)-B(119)
! JVS(101) = Jac_FULL(53,54)
  JVS(101) = B(11)
! JVS(102) = Jac_FULL(54,53)
  JVS(102) = B(23)
! JVS(103) = Jac_FULL(54,54)
  JVS(103) = -B(11)
! JVS(104) = Jac_FULL(55,55)
  JVS(104) = -B(24)-B(68)-B(120)
! JVS(105) = Jac_FULL(55,56)
  JVS(105) = B(12)
! JVS(106) = Jac_FULL(56,55)
  JVS(106) = B(24)
! JVS(107) = Jac_FULL(56,56)
  JVS(107) = -B(12)
! JVS(108) = Jac_FULL(57,57)
  JVS(108) = -B(9)
! JVS(109) = Jac_FULL(57,75)
  JVS(109) = B(21)
! JVS(110) = Jac_FULL(58,58)
  JVS(110) = -B(10)
! JVS(111) = Jac_FULL(58,68)
  JVS(111) = B(22)
! JVS(112) = Jac_FULL(59,59)
  JVS(112) = -B(1)
! JVS(113) = Jac_FULL(59,82)
  JVS(113) = B(13)
! JVS(114) = Jac_FULL(60,46)
  JVS(114) = B(2)
! JVS(115) = Jac_FULL(60,60)
  JVS(115) = -B(14)-B(35)-B(109)
! JVS(116) = Jac_FULL(60,61)
  JVS(116) = B(48)
! JVS(117) = Jac_FULL(60,79)
  JVS(117) = B(49)
! JVS(118) = Jac_FULL(61,1)
  JVS(118) = B(25)
! JVS(119) = Jac_FULL(61,60)
  JVS(119) = B(35)
! JVS(120) = Jac_FULL(61,61)
  JVS(120) = -B(48)-B(122)
! JVS(121) = Jac_FULL(61,79)
  JVS(121) = -B(49)
! JVS(122) = Jac_FULL(62,48)
  JVS(122) = B(4)
! JVS(123) = Jac_FULL(62,62)
  JVS(123) = -B(16)-B(38)-B(111)
! JVS(124) = Jac_FULL(62,63)
  JVS(124) = B(54)
! JVS(125) = Jac_FULL(62,64)
  JVS(125) = B(55)
! JVS(126) = Jac_FULL(63,2)
  JVS(126) = B(26)
! JVS(127) = Jac_FULL(63,62)
  JVS(127) = B(38)
! JVS(128) = Jac_FULL(63,63)
  JVS(128) = -B(54)-B(123)
! JVS(129) = Jac_FULL(63,64)
  JVS(129) = -B(55)
! JVS(130) = Jac_FULL(64,62)
  JVS(130) = B(38)
! JVS(131) = Jac_FULL(64,63)
  JVS(131) = -B(54)
! JVS(132) = Jac_FULL(64,64)
  JVS(132) = -B(55)-B(60)-B(136)
! JVS(133) = Jac_FULL(64,79)
  JVS(133) = -B(61)
! JVS(134) = Jac_FULL(65,51)
  JVS(134) = B(7)
! JVS(135) = Jac_FULL(65,65)
  JVS(135) = -B(19)-B(39)-B(114)
! JVS(136) = Jac_FULL(65,66)
  JVS(136) = B(56)
! JVS(137) = Jac_FULL(65,79)
  JVS(137) = B(57)
! JVS(138) = Jac_FULL(66,65)
  JVS(138) = B(39)
! JVS(139) = Jac_FULL(66,66)
  JVS(139) = -B(56)-B(141)
! JVS(140) = Jac_FULL(66,79)
  JVS(140) = -B(57)
! JVS(141) = Jac_FULL(67,67)
  JVS(141) = -B(42)-B(117)
! JVS(142) = Jac_FULL(67,73)
  JVS(142) = B(62)
! JVS(143) = Jac_FULL(67,79)
  JVS(143) = B(63)
! JVS(144) = Jac_FULL(68,58)
  JVS(144) = B(10)
! JVS(145) = Jac_FULL(68,68)
  JVS(145) = -B(22)-B(40)-B(118)
! JVS(146) = Jac_FULL(68,69)
  JVS(146) = B(58)
! JVS(147) = Jac_FULL(68,79)
  JVS(147) = B(59)
! JVS(148) = Jac_FULL(69,3)
  JVS(148) = B(27)
! JVS(149) = Jac_FULL(69,68)
  JVS(149) = B(40)
! JVS(150) = Jac_FULL(69,69)
  JVS(150) = -B(58)-B(124)
! JVS(151) = Jac_FULL(69,79)
  JVS(151) = -B(59)
! JVS(152) = Jac_FULL(70,70)
  JVS(152) = -B(52)-B(140)
! JVS(153) = Jac_FULL(70,72)
  JVS(153) = B(37)
! JVS(154) = Jac_FULL(70,79)
  JVS(154) = -B(53)
! JVS(155) = Jac_FULL(71,47)
  JVS(155) = B(3)
! JVS(156) = Jac_FULL(71,71)
  JVS(156) = -B(15)-B(36)-B(110)
! JVS(157) = Jac_FULL(71,72)
  JVS(157) = B(50)
! JVS(158) = Jac_FULL(71,79)
  JVS(158) = B(51)
! JVS(159) = Jac_FULL(72,70)
  JVS(159) = B(52)
! JVS(160) = Jac_FULL(72,71)
  JVS(160) = B(36)
! JVS(161) = Jac_FULL(72,72)
  JVS(161) = -B(37)-B(50)-B(139)
! JVS(162) = Jac_FULL(72,79)
  JVS(162) = -B(51)+B(53)
! JVS(163) = Jac_FULL(73,67)
  JVS(163) = B(42)
! JVS(164) = Jac_FULL(73,73)
  JVS(164) = -B(43)-B(62)-B(142)
! JVS(165) = Jac_FULL(73,76)
  JVS(165) = B(64)
! JVS(166) = Jac_FULL(73,79)
  JVS(166) = -B(63)+B(65)
! JVS(167) = Jac_FULL(74,49)
  JVS(167) = B(5)
! JVS(168) = Jac_FULL(74,74)
  JVS(168) = -B(17)-B(70)-B(112)
! JVS(169) = Jac_FULL(74,79)
  JVS(169) = -B(71)
! JVS(170) = Jac_FULL(74,81)
  JVS(170) = -B(72)
! JVS(171) = Jac_FULL(75,57)
  JVS(171) = B(9)
! JVS(172) = Jac_FULL(75,75)
  JVS(172) = -B(21)-B(103)-B(106)-B(116)
! JVS(173) = Jac_FULL(75,79)
  JVS(173) = -B(104)
! JVS(174) = Jac_FULL(75,81)
  JVS(174) = -B(105)-B(107)
! JVS(175) = Jac_FULL(76,5)
  JVS(175) = B(29)
! JVS(176) = Jac_FULL(76,12)
  JVS(176) = B(85)+B(87)+B(89)+B(91)+B(94)+B(97)
! JVS(177) = Jac_FULL(76,13)
  JVS(177) = B(79)+B(81)+B(83)+B(92)+B(95)+B(98)
! JVS(178) = Jac_FULL(76,73)
  JVS(178) = B(43)
! JVS(179) = Jac_FULL(76,74)
  JVS(179) = B(70)
! JVS(180) = Jac_FULL(76,75)
  JVS(180) = B(103)+B(106)
! JVS(181) = Jac_FULL(76,76)
  JVS(181) = -B(64)-B(121)
! JVS(182) = Jac_FULL(76,77)
  JVS(182) = B(100)
! JVS(183) = Jac_FULL(76,78)
  JVS(183) = B(77)+B(84)+B(90)+B(99)
! JVS(184) = Jac_FULL(76,79)
  JVS(184) = -B(65)+B(71)+B(101)+B(104)
! JVS(185) = Jac_FULL(76,80)
  JVS(185) = B(73)+B(75)+B(78)
! JVS(186) = Jac_FULL(76,81)
  JVS(186) = B(72)+B(76)+B(82)+B(88)+B(96)+B(102)+B(105)+B(107)
! JVS(187) = Jac_FULL(76,82)
  JVS(187) = B(74)+B(80)+B(86)+B(93)
! JVS(188) = Jac_FULL(77,52)
  JVS(188) = B(8)
! JVS(189) = Jac_FULL(77,77)
  JVS(189) = -B(20)-B(100)-B(115)
! JVS(190) = Jac_FULL(77,79)
  JVS(190) = -B(101)
! JVS(191) = Jac_FULL(77,81)
  JVS(191) = -B(102)
! JVS(192) = Jac_FULL(78,12)
  JVS(192) = -B(89)-B(97)
! JVS(193) = Jac_FULL(78,13)
  JVS(193) = -B(83)-B(98)
! JVS(194) = Jac_FULL(78,78)
  JVS(194) = -B(46)-B(77)-B(84)-B(90)-B(99)-B(138)
! JVS(195) = Jac_FULL(78,79)
  JVS(195) = -B(47)
! JVS(196) = Jac_FULL(78,80)
  JVS(196) = -B(78)
! JVS(197) = Jac_FULL(78,81)
  JVS(197) = B(34)
! JVS(198) = Jac_FULL(79,12)
  JVS(198) = 2*B(85)+B(87)+2*B(91)+B(94)
! JVS(199) = Jac_FULL(79,13)
  JVS(199) = 2*B(79)+B(81)+2*B(92)+B(95)
! JVS(200) = Jac_FULL(79,60)
  JVS(200) = B(35)
! JVS(201) = Jac_FULL(79,61)
  JVS(201) = -B(48)
! JVS(202) = Jac_FULL(79,64)
  JVS(202) = -B(60)
! JVS(203) = Jac_FULL(79,65)
  JVS(203) = B(39)
! JVS(204) = Jac_FULL(79,66)
  JVS(204) = -B(56)
! JVS(205) = Jac_FULL(79,67)
  JVS(205) = B(42)
! JVS(206) = Jac_FULL(79,68)
  JVS(206) = B(40)
! JVS(207) = Jac_FULL(79,69)
  JVS(207) = -B(58)
! JVS(208) = Jac_FULL(79,70)
  JVS(208) = -B(52)
! JVS(209) = Jac_FULL(79,71)
  JVS(209) = B(36)
! JVS(210) = Jac_FULL(79,72)
  JVS(210) = B(37)-B(50)
! JVS(211) = Jac_FULL(79,73)
  JVS(211) = B(43)-B(62)
! JVS(212) = Jac_FULL(79,74)
  JVS(212) = B(70)
! JVS(213) = Jac_FULL(79,75)
  JVS(213) = B(103)+B(106)
! JVS(214) = Jac_FULL(79,76)
  JVS(214) = -B(64)
! JVS(215) = Jac_FULL(79,77)
  JVS(215) = B(100)
! JVS(216) = Jac_FULL(79,78)
  JVS(216) = -B(46)
! JVS(217) = Jac_FULL(79,79)
  JVS(217) = -B(44)-B(47)-B(49)-B(51)-B(53)-B(57)-B(59)-B(61)-B(63)-B(65)+B(71)+B(101)+B(104)-B(135)
! JVS(218) = Jac_FULL(79,80)
  JVS(218) = 2*B(73)+B(75)
! JVS(219) = Jac_FULL(79,81)
  JVS(219) = B(34)-B(45)+B(72)+B(76)+B(82)+B(88)+B(96)+B(102)+B(105)+B(107)
! JVS(220) = Jac_FULL(79,82)
  JVS(220) = B(33)+2*B(74)+2*B(80)+2*B(86)+2*B(93)
! JVS(221) = Jac_FULL(80,50)
  JVS(221) = B(6)
! JVS(222) = Jac_FULL(80,78)
  JVS(222) = -B(77)
! JVS(223) = Jac_FULL(80,79)
  JVS(223) = 0
! JVS(224) = Jac_FULL(80,80)
  JVS(224) = -B(18)-B(73)-B(75)-B(78)-B(113)
! JVS(225) = Jac_FULL(80,81)
  JVS(225) = -B(76)
! JVS(226) = Jac_FULL(80,82)
  JVS(226) = -B(74)
! JVS(227) = Jac_FULL(81,12)
  JVS(227) = -B(87)-B(94)
! JVS(228) = Jac_FULL(81,13)
  JVS(228) = -B(81)-B(95)
! JVS(229) = Jac_FULL(81,74)
  JVS(229) = -B(70)
! JVS(230) = Jac_FULL(81,75)
  JVS(230) = -B(103)-B(106)
! JVS(231) = Jac_FULL(81,77)
  JVS(231) = -B(100)
! JVS(232) = Jac_FULL(81,78)
  JVS(232) = B(46)
! JVS(233) = Jac_FULL(81,79)
  JVS(233) = -B(44)+B(47)-B(71)-B(101)-B(104)
! JVS(234) = Jac_FULL(81,80)
  JVS(234) = -B(75)
! JVS(235) = Jac_FULL(81,81)
  JVS(235) = -B(34)-B(45)-B(72)-B(76)-B(82)-B(88)-B(96)-B(102)-B(105)-B(107)-B(137)
! JVS(236) = Jac_FULL(81,82)
  JVS(236) = B(33)
! JVS(237) = Jac_FULL(82,12)
  JVS(237) = -B(85)-B(91)
! JVS(238) = Jac_FULL(82,13)
  JVS(238) = -B(79)-B(92)
! JVS(239) = Jac_FULL(82,59)
  JVS(239) = B(1)
! JVS(240) = Jac_FULL(82,79)
  JVS(240) = B(44)
! JVS(241) = Jac_FULL(82,80)
  JVS(241) = -B(73)
! JVS(242) = Jac_FULL(82,81)
  JVS(242) = B(45)
! JVS(243) = Jac_FULL(82,82)
  JVS(243) = -B(13)-B(33)-B(74)-B(80)-B(86)-B(93)-B(108)
      
END SUBROUTINE Jac_SP

! End of Jac_SP function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Jac_SP_Vec - function for sparse multiplication: sparse Jacobian times vector
!   Arguments :
!      JVS       - sparse Jacobian of variables
!      UV        - User vector for variables
!      JUV       - Jacobian times user vector
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Jac_SP_Vec ( JVS, UV, JUV )

! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)
! UV - User vector for variables
  REAL(kind=dp) :: UV(NVAR)
! JUV - Jacobian times user vector
  REAL(kind=dp) :: JUV(NVAR)

  JUV(1) = JVS(1)*UV(1)
  JUV(2) = JVS(2)*UV(2)
  JUV(3) = JVS(3)*UV(3)
  JUV(4) = JVS(4)*UV(4)
  JUV(5) = JVS(5)*UV(5)
  JUV(6) = JVS(6)*UV(6)
  JUV(7) = JVS(7)*UV(7)
  JUV(8) = JVS(8)*UV(8)
  JUV(9) = JVS(9)*UV(8)+JVS(10)*UV(9)
  JUV(10) = JVS(11)*UV(6)+JVS(12)*UV(10)
  JUV(11) = JVS(13)*UV(7)+JVS(14)*UV(11)
  JUV(12) = JVS(15)*UV(12)
  JUV(13) = JVS(16)*UV(13)
  JUV(14) = JVS(17)*UV(14)
  JUV(15) = JVS(18)*UV(15)
  JUV(16) = JVS(19)*UV(16)
  JUV(17) = JVS(20)*UV(4)+JVS(21)*UV(17)
  JUV(18) = JVS(22)*UV(18)+JVS(23)*UV(78)+JVS(24)*UV(81)+JVS(25)*UV(82)
  JUV(19) = JVS(26)*UV(19)+JVS(27)*UV(60)
  JUV(20) = JVS(28)*UV(20)+JVS(29)*UV(70)+JVS(30)*UV(71)+JVS(31)*UV(72)
  JUV(21) = JVS(32)*UV(21)+JVS(33)*UV(62)
  JUV(22) = JVS(34)*UV(22)+JVS(35)*UV(74)
  JUV(23) = JVS(36)*UV(23)+JVS(37)*UV(80)
  JUV(24) = JVS(38)*UV(24)+JVS(39)*UV(65)+JVS(40)*UV(66)
  JUV(25) = JVS(41)*UV(25)+JVS(42)*UV(77)
  JUV(26) = JVS(43)*UV(26)+JVS(44)*UV(75)
  JUV(27) = JVS(45)*UV(27)+JVS(46)*UV(67)+JVS(47)*UV(73)+JVS(48)*UV(76)
  JUV(28) = JVS(49)*UV(28)+JVS(50)*UV(68)
  JUV(29) = JVS(51)*UV(29)+JVS(52)*UV(53)
  JUV(30) = JVS(53)*UV(30)+JVS(54)*UV(55)
  JUV(31) = JVS(55)*UV(31)+JVS(56)*UV(61)
  JUV(32) = JVS(57)*UV(32)+JVS(58)*UV(63)
  JUV(33) = JVS(59)*UV(33)+JVS(60)*UV(69)
  JUV(34) = JVS(61)*UV(9)+JVS(62)*UV(34)
  JUV(35) = JVS(63)*UV(12)+JVS(64)*UV(35)
  JUV(36) = JVS(65)*UV(13)+JVS(66)*UV(36)
  JUV(37) = JVS(67)*UV(10)+JVS(68)*UV(37)
  JUV(38) = JVS(69)*UV(38)+JVS(70)*UV(39)
  JUV(39) = JVS(71)*UV(39)+JVS(72)*UV(53)+JVS(73)*UV(55)
  JUV(40) = JVS(74)*UV(11)+JVS(75)*UV(40)
  JUV(41) = JVS(76)*UV(41)+JVS(77)*UV(79)
  JUV(42) = JVS(78)*UV(14)+JVS(79)*UV(42)
  JUV(43) = JVS(80)*UV(15)+JVS(81)*UV(43)
  JUV(44) = JVS(82)*UV(16)+JVS(83)*UV(44)
  JUV(45) = JVS(84)*UV(17)+JVS(85)*UV(45)
  JUV(46) = JVS(86)*UV(46)+JVS(87)*UV(60)
  JUV(47) = JVS(88)*UV(47)+JVS(89)*UV(71)
  JUV(48) = JVS(90)*UV(48)+JVS(91)*UV(62)
  JUV(49) = JVS(92)*UV(49)+JVS(93)*UV(74)
  JUV(50) = JVS(94)*UV(50)+JVS(95)*UV(80)
  JUV(51) = JVS(96)*UV(51)+JVS(97)*UV(65)
  JUV(52) = JVS(98)*UV(52)+JVS(99)*UV(77)
  JUV(53) = JVS(100)*UV(53)+JVS(101)*UV(54)
  JUV(54) = JVS(102)*UV(53)+JVS(103)*UV(54)
  JUV(55) = JVS(104)*UV(55)+JVS(105)*UV(56)
  JUV(56) = JVS(106)*UV(55)+JVS(107)*UV(56)
  JUV(57) = JVS(108)*UV(57)+JVS(109)*UV(75)
  JUV(58) = JVS(110)*UV(58)+JVS(111)*UV(68)
  JUV(59) = JVS(112)*UV(59)+JVS(113)*UV(82)
  JUV(60) = JVS(114)*UV(46)+JVS(115)*UV(60)+JVS(116)*UV(61)+JVS(117)*UV(79)
  JUV(61) = JVS(118)*UV(1)+JVS(119)*UV(60)+JVS(120)*UV(61)+JVS(121)*UV(79)
  JUV(62) = JVS(122)*UV(48)+JVS(123)*UV(62)+JVS(124)*UV(63)+JVS(125)*UV(64)
  JUV(63) = JVS(126)*UV(2)+JVS(127)*UV(62)+JVS(128)*UV(63)+JVS(129)*UV(64)
  JUV(64) = JVS(130)*UV(62)+JVS(131)*UV(63)+JVS(132)*UV(64)+JVS(133)*UV(79)
  JUV(65) = JVS(134)*UV(51)+JVS(135)*UV(65)+JVS(136)*UV(66)+JVS(137)*UV(79)
  JUV(66) = JVS(138)*UV(65)+JVS(139)*UV(66)+JVS(140)*UV(79)
  JUV(67) = JVS(141)*UV(67)+JVS(142)*UV(73)+JVS(143)*UV(79)
  JUV(68) = JVS(144)*UV(58)+JVS(145)*UV(68)+JVS(146)*UV(69)+JVS(147)*UV(79)
  JUV(69) = JVS(148)*UV(3)+JVS(149)*UV(68)+JVS(150)*UV(69)+JVS(151)*UV(79)
  JUV(70) = JVS(152)*UV(70)+JVS(153)*UV(72)+JVS(154)*UV(79)
  JUV(71) = JVS(155)*UV(47)+JVS(156)*UV(71)+JVS(157)*UV(72)+JVS(158)*UV(79)
  JUV(72) = JVS(159)*UV(70)+JVS(160)*UV(71)+JVS(161)*UV(72)+JVS(162)*UV(79)
  JUV(73) = JVS(163)*UV(67)+JVS(164)*UV(73)+JVS(165)*UV(76)+JVS(166)*UV(79)
  JUV(74) = JVS(167)*UV(49)+JVS(168)*UV(74)+JVS(169)*UV(79)+JVS(170)*UV(81)
  JUV(75) = JVS(171)*UV(57)+JVS(172)*UV(75)+JVS(173)*UV(79)+JVS(174)*UV(81)
  JUV(76) = JVS(175)*UV(5)+JVS(176)*UV(12)+JVS(177)*UV(13)+JVS(178)*UV(73)+JVS(179)*UV(74)+JVS(180)*UV(75)+JVS(181)&
              &*UV(76)+JVS(182)*UV(77)+JVS(183)*UV(78)+JVS(184)*UV(79)+JVS(185)*UV(80)+JVS(186)*UV(81)+JVS(187)*UV(82)
  JUV(77) = JVS(188)*UV(52)+JVS(189)*UV(77)+JVS(190)*UV(79)+JVS(191)*UV(81)
  JUV(78) = JVS(192)*UV(12)+JVS(193)*UV(13)+JVS(194)*UV(78)+JVS(195)*UV(79)+JVS(196)*UV(80)+JVS(197)*UV(81)
  JUV(79) = JVS(198)*UV(12)+JVS(199)*UV(13)+JVS(200)*UV(60)+JVS(201)*UV(61)+JVS(202)*UV(64)+JVS(203)*UV(65)+JVS(204)&
              &*UV(66)+JVS(205)*UV(67)+JVS(206)*UV(68)+JVS(207)*UV(69)+JVS(208)*UV(70)+JVS(209)*UV(71)+JVS(210)*UV(72)&
              &+JVS(211)*UV(73)+JVS(212)*UV(74)+JVS(213)*UV(75)+JVS(214)*UV(76)+JVS(215)*UV(77)+JVS(216)*UV(78)+JVS(217)&
              &*UV(79)+JVS(218)*UV(80)+JVS(219)*UV(81)+JVS(220)*UV(82)
  JUV(80) = JVS(221)*UV(50)+JVS(222)*UV(78)+JVS(224)*UV(80)+JVS(225)*UV(81)+JVS(226)*UV(82)
  JUV(81) = JVS(227)*UV(12)+JVS(228)*UV(13)+JVS(229)*UV(74)+JVS(230)*UV(75)+JVS(231)*UV(77)+JVS(232)*UV(78)+JVS(233)&
              &*UV(79)+JVS(234)*UV(80)+JVS(235)*UV(81)+JVS(236)*UV(82)
  JUV(82) = JVS(237)*UV(12)+JVS(238)*UV(13)+JVS(239)*UV(59)+JVS(240)*UV(79)+JVS(241)*UV(80)+JVS(242)*UV(81)+JVS(243)&
              &*UV(82)
      
END SUBROUTINE Jac_SP_Vec

! End of Jac_SP_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! JacTR_SP_Vec - sparse multiplication: sparse Jacobian transposed times vector
!   Arguments :
!      JVS       - sparse Jacobian of variables
!      UV        - User vector for variables
!      JTUV      - Jacobian transposed times user vector
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE JacTR_SP_Vec ( JVS, UV, JTUV )

! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)
! UV - User vector for variables
  REAL(kind=dp) :: UV(NVAR)
! JTUV - Jacobian transposed times user vector
  REAL(kind=dp) :: JTUV(NVAR)

  JTUV(1) = JVS(1)*UV(1)+JVS(118)*UV(61)
  JTUV(2) = JVS(2)*UV(2)+JVS(126)*UV(63)
  JTUV(3) = JVS(3)*UV(3)+JVS(148)*UV(69)
  JTUV(4) = JVS(4)*UV(4)+JVS(20)*UV(17)
  JTUV(5) = JVS(5)*UV(5)+JVS(175)*UV(76)
  JTUV(6) = JVS(6)*UV(6)+JVS(11)*UV(10)
  JTUV(7) = JVS(7)*UV(7)+JVS(13)*UV(11)
  JTUV(8) = JVS(8)*UV(8)+JVS(9)*UV(9)
  JTUV(9) = JVS(10)*UV(9)+JVS(61)*UV(34)
  JTUV(10) = JVS(12)*UV(10)+JVS(67)*UV(37)
  JTUV(11) = JVS(14)*UV(11)+JVS(74)*UV(40)
  JTUV(12) = JVS(15)*UV(12)+JVS(63)*UV(35)+JVS(176)*UV(76)+JVS(192)*UV(78)+JVS(198)*UV(79)+JVS(227)*UV(81)+JVS(237)&
               &*UV(82)
  JTUV(13) = JVS(16)*UV(13)+JVS(65)*UV(36)+JVS(177)*UV(76)+JVS(193)*UV(78)+JVS(199)*UV(79)+JVS(228)*UV(81)+JVS(238)&
               &*UV(82)
  JTUV(14) = JVS(17)*UV(14)+JVS(78)*UV(42)
  JTUV(15) = JVS(18)*UV(15)+JVS(80)*UV(43)
  JTUV(16) = JVS(19)*UV(16)+JVS(82)*UV(44)
  JTUV(17) = JVS(21)*UV(17)+JVS(84)*UV(45)
  JTUV(18) = JVS(22)*UV(18)
  JTUV(19) = JVS(26)*UV(19)
  JTUV(20) = JVS(28)*UV(20)
  JTUV(21) = JVS(32)*UV(21)
  JTUV(22) = JVS(34)*UV(22)
  JTUV(23) = JVS(36)*UV(23)
  JTUV(24) = JVS(38)*UV(24)
  JTUV(25) = JVS(41)*UV(25)
  JTUV(26) = JVS(43)*UV(26)
  JTUV(27) = JVS(45)*UV(27)
  JTUV(28) = JVS(49)*UV(28)
  JTUV(29) = JVS(51)*UV(29)
  JTUV(30) = JVS(53)*UV(30)
  JTUV(31) = JVS(55)*UV(31)
  JTUV(32) = JVS(57)*UV(32)
  JTUV(33) = JVS(59)*UV(33)
  JTUV(34) = JVS(62)*UV(34)
  JTUV(35) = JVS(64)*UV(35)
  JTUV(36) = JVS(66)*UV(36)
  JTUV(37) = JVS(68)*UV(37)
  JTUV(38) = JVS(69)*UV(38)
  JTUV(39) = JVS(70)*UV(38)+JVS(71)*UV(39)
  JTUV(40) = JVS(75)*UV(40)
  JTUV(41) = JVS(76)*UV(41)
  JTUV(42) = JVS(79)*UV(42)
  JTUV(43) = JVS(81)*UV(43)
  JTUV(44) = JVS(83)*UV(44)
  JTUV(45) = JVS(85)*UV(45)
  JTUV(46) = JVS(86)*UV(46)+JVS(114)*UV(60)
  JTUV(47) = JVS(88)*UV(47)+JVS(155)*UV(71)
  JTUV(48) = JVS(90)*UV(48)+JVS(122)*UV(62)
  JTUV(49) = JVS(92)*UV(49)+JVS(167)*UV(74)
  JTUV(50) = JVS(94)*UV(50)+JVS(221)*UV(80)
  JTUV(51) = JVS(96)*UV(51)+JVS(134)*UV(65)
  JTUV(52) = JVS(98)*UV(52)+JVS(188)*UV(77)
  JTUV(53) = JVS(52)*UV(29)+JVS(72)*UV(39)+JVS(100)*UV(53)+JVS(102)*UV(54)
  JTUV(54) = JVS(101)*UV(53)+JVS(103)*UV(54)
  JTUV(55) = JVS(54)*UV(30)+JVS(73)*UV(39)+JVS(104)*UV(55)+JVS(106)*UV(56)
  JTUV(56) = JVS(105)*UV(55)+JVS(107)*UV(56)
  JTUV(57) = JVS(108)*UV(57)+JVS(171)*UV(75)
  JTUV(58) = JVS(110)*UV(58)+JVS(144)*UV(68)
  JTUV(59) = JVS(112)*UV(59)+JVS(239)*UV(82)
  JTUV(60) = JVS(27)*UV(19)+JVS(87)*UV(46)+JVS(115)*UV(60)+JVS(119)*UV(61)+JVS(200)*UV(79)
  JTUV(61) = JVS(56)*UV(31)+JVS(116)*UV(60)+JVS(120)*UV(61)+JVS(201)*UV(79)
  JTUV(62) = JVS(33)*UV(21)+JVS(91)*UV(48)+JVS(123)*UV(62)+JVS(127)*UV(63)+JVS(130)*UV(64)
  JTUV(63) = JVS(58)*UV(32)+JVS(124)*UV(62)+JVS(128)*UV(63)+JVS(131)*UV(64)
  JTUV(64) = JVS(125)*UV(62)+JVS(129)*UV(63)+JVS(132)*UV(64)+JVS(202)*UV(79)
  JTUV(65) = JVS(39)*UV(24)+JVS(97)*UV(51)+JVS(135)*UV(65)+JVS(138)*UV(66)+JVS(203)*UV(79)
  JTUV(66) = JVS(40)*UV(24)+JVS(136)*UV(65)+JVS(139)*UV(66)+JVS(204)*UV(79)
  JTUV(67) = JVS(46)*UV(27)+JVS(141)*UV(67)+JVS(163)*UV(73)+JVS(205)*UV(79)
  JTUV(68) = JVS(50)*UV(28)+JVS(111)*UV(58)+JVS(145)*UV(68)+JVS(149)*UV(69)+JVS(206)*UV(79)
  JTUV(69) = JVS(60)*UV(33)+JVS(146)*UV(68)+JVS(150)*UV(69)+JVS(207)*UV(79)
  JTUV(70) = JVS(29)*UV(20)+JVS(152)*UV(70)+JVS(159)*UV(72)+JVS(208)*UV(79)
  JTUV(71) = JVS(30)*UV(20)+JVS(89)*UV(47)+JVS(156)*UV(71)+JVS(160)*UV(72)+JVS(209)*UV(79)
  JTUV(72) = JVS(31)*UV(20)+JVS(153)*UV(70)+JVS(157)*UV(71)+JVS(161)*UV(72)+JVS(210)*UV(79)
  JTUV(73) = JVS(47)*UV(27)+JVS(142)*UV(67)+JVS(164)*UV(73)+JVS(178)*UV(76)+JVS(211)*UV(79)
  JTUV(74) = JVS(35)*UV(22)+JVS(93)*UV(49)+JVS(168)*UV(74)+JVS(179)*UV(76)+JVS(212)*UV(79)+JVS(229)*UV(81)
  JTUV(75) = JVS(44)*UV(26)+JVS(109)*UV(57)+JVS(172)*UV(75)+JVS(180)*UV(76)+JVS(213)*UV(79)+JVS(230)*UV(81)
  JTUV(76) = JVS(48)*UV(27)+JVS(165)*UV(73)+JVS(181)*UV(76)+JVS(214)*UV(79)
  JTUV(77) = JVS(42)*UV(25)+JVS(99)*UV(52)+JVS(182)*UV(76)+JVS(189)*UV(77)+JVS(215)*UV(79)+JVS(231)*UV(81)
  JTUV(78) = JVS(23)*UV(18)+JVS(183)*UV(76)+JVS(194)*UV(78)+JVS(216)*UV(79)+JVS(222)*UV(80)+JVS(232)*UV(81)
  JTUV(79) = JVS(77)*UV(41)+JVS(117)*UV(60)+JVS(121)*UV(61)+JVS(133)*UV(64)+JVS(137)*UV(65)+JVS(140)*UV(66)+JVS(143)&
               &*UV(67)+JVS(147)*UV(68)+JVS(151)*UV(69)+JVS(154)*UV(70)+JVS(158)*UV(71)+JVS(162)*UV(72)+JVS(166)*UV(73)&
               &+JVS(169)*UV(74)+JVS(173)*UV(75)+JVS(184)*UV(76)+JVS(190)*UV(77)+JVS(195)*UV(78)+JVS(217)*UV(79)+JVS(233)&
               &*UV(81)+JVS(240)*UV(82)
  JTUV(80) = JVS(37)*UV(23)+JVS(95)*UV(50)+JVS(185)*UV(76)+JVS(196)*UV(78)+JVS(218)*UV(79)+JVS(224)*UV(80)+JVS(234)&
               &*UV(81)+JVS(241)*UV(82)
  JTUV(81) = JVS(24)*UV(18)+JVS(170)*UV(74)+JVS(174)*UV(75)+JVS(186)*UV(76)+JVS(191)*UV(77)+JVS(197)*UV(78)+JVS(219)&
               &*UV(79)+JVS(225)*UV(80)+JVS(235)*UV(81)+JVS(242)*UV(82)
  JTUV(82) = JVS(25)*UV(18)+JVS(113)*UV(59)+JVS(187)*UV(76)+JVS(220)*UV(79)+JVS(226)*UV(80)+JVS(236)*UV(81)+JVS(243)&
               &*UV(82)
      
END SUBROUTINE JacTR_SP_Vec

! End of JacTR_SP_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE aqchem_Jacobian

