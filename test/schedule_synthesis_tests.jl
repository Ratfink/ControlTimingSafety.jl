using LinearAlgebra
using RealTimeScheduling

@testset "Schedule Synthesis" begin
     sys = let
          r_1 = 100000
          r_2 = 500000
          r_3 = 200000
          c_1 = 0.000002
          c_2 = 0.000010
          A = [-1/c_1 * (1/r_1 + 1/r_2)  1/(r_2*c_1)
               1/(r_2*c_2)               -1/c_2 * (1/r_2 + 1/r_3)]
          B = [1/(r_1*c_1)
               1/(r_3*c_2)]
          #C = [1  -1]
          C = [1  0
               0  1]
          D = 0
          ss(A, B, C, D)
     end
     ctrl_delay = 0.1
     sysd = c2d(sys, ctrl_delay)
     sysd_delay = let
          A = [sysd.A sysd.B; 0 0 0]
          B = [0; 0; 1]
          C = [sysd.C zeros(2, 1)]
          D = sysd.D
          ss(A, B, C, D, sysd.Ts)
     end
     K = let
          state_cost = 2
          control_cost = 1
          Q = I * state_cost
          R = I * control_cost
          lqr(sysd_delay, Q, R)
     end

     @test synthesize_constraints(sysd, K, repeat([10.0], 3), 5.5, 6, 10, 10) == [
          RealTimeScheduling.MeetAny(1, 2),
          RealTimeScheduling.MeetAny(2, 3),
          RealTimeScheduling.MeetAny(3, 4),
          RealTimeScheduling.MeetAny(4, 5),
          RealTimeScheduling.MeetAny(5, 6)
     ]

     @test schedule_xghtc([
          RealTimeScheduling.MeetAny(1, 3),
          RealTimeScheduling.MeetAny(2, 3),
          RealTimeScheduling.MeetAny(1, 4),
          RealTimeScheduling.MeetAny(1, 4),
          RealTimeScheduling.MeetAny(1, 2)
     ], slotsize=2) == [
          # 0  0  1  0  0  1  0  0  1  0  0  1  0  0  1  0
          # 1  1  0  1  1  0  1  1  0  1  1  0  1  1  0  1
          # 1  0  0  0  1  0  0  0  1  0  0  0  1  0  0  0
          # 0  0  1  0  0  0  1  0  0  0  1  0  0  0  1  0
          # 0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1
          0  1  0  0  1  0  0  1  0  0  1  0
          1  0  1  1  0  1  1  0  1  1  0  1
          1  0  0  0  1  0  0  0  1  0  0  0
          0  0  1  0  0  0  1  0  0  0  1  0
          0  1  0  1  0  1  0  1  0  1  0  1
     ]
end
