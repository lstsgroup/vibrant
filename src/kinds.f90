!
!   Copyright 2025 Ekin E. Winogradow, Johannes Scheffler, Moritz Leucke, Dorothea Golze
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.
!

!> @brief Module containing the kind definitions for the real numbers
MODULE kinds

    IMPLICIT NONE

    !> Double precision kind
    INTEGER, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND(15)

    !> Single precision kind
    INTEGER, PARAMETER, PUBLIC :: sp = SELECTED_REAL_KIND(6)

    !> default string length
    INTEGER, PARAMETER, PUBLIC :: str_len = 100

END MODULE kinds
