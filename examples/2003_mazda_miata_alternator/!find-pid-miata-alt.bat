rem Alternator PID for 2003 Mazda Miata:

@rem 54=debugFloatField1=altpid.output
@rem 10=VBatt
@rem 1=rpm
@rem 22=IAC_position
@rem 12=advance angle

@rem Usage: PID_FROM_MSL <file.msl> <startTime> <endTime> <inputColumnIdx> <outputColumnIdx> [<targetValue>]
..\..\pid_from_msl.exe "miata_alt_log.msl" 66.737 71.662 54 10 13.8
