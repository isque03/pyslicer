M109 S200.0 ; Heat up to 200.0C
G90       ; Use absolute coordinates
G21       ; Set units to millimeters
; pyslicer planning: outer=70.0mm/s inner=70.0mm/s infill=70.0mm/s accel=1000mm/s^2 jerk=20.0mm/s corner=5.0mm/s min_angle=20deg
M106 S0   ; Fan Off
G28       ; Home all axes
G92 E0    ; Zero extruder
M82       ; Use absolute distances for extrusion
G92 E0    ; Zero extruder
G1 F200 E8  ; prime extruder
G92 E0    ; Zero extruder
M117 Printing.

;; New Layer Z: 0.2 
;; Perimeter 0
G1 F2400.000000 Z0.200000
;; Contour 0 Area: 100.0
G1 F4200.000000
;; travel move
G1 X0.000000 Y0.000000
G1 F300.000000 X10.000000 Y0.000000 E0.332601
G1 X10.000000 Y10.000000 E0.665203
G1 X0.000000 Y10.000000 E0.997804
G1 X0.000000 Y0.000000 E1.330405
;; Infill
;;infill travel move. distance: 0.000000 
G1 F8000.000000 X1.000000 Y1.000000
G1 F4200.000000
G1 F2683.281573 X2.000000 Y1.000000 E1.363666
G1 X3.000000 Y1.000000 E1.396926
;; retract 
;;infill travel move. distance: 6.403124 
G1 F8000.000000 X8.000000 Y5.000000
G1 F4200.000000
G1 X9.000000 Y5.000000 E1.430186
M107    ; Fan off
M104 S0 ; Heat off 200.0C
M140 S0 ; Bed heat off
G91     ; Relative
G1 E-1 F400 ; Retract
G1 Z+1.0 E-5 X-20 Y-20 F9000
G28 X0 Y0
M84     ; Motors off
G90     ; Absolute
