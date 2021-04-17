        select case (this%periodic)
        case (.TRUE.)
                RHS(:,:,         1) = ao2Penta * ( f(:,:,         1) )                     &
                                    + bo2Penta * ( f(:,:,         2) + f(:,:,    this%n) ) &
                                    + co2Penta * ( f(:,:,         3) + f(:,:,  this%n-1) ) &
                                    + do2Penta * ( f(:,:,         4) + f(:,:,  this%n-2) ) &
                                    + eo2Penta * ( f(:,:,         5) + f(:,:,  this%n-3) )
                RHS(:,:,         2) = ao2Penta * ( f(:,:,         2) )                     &
                                    + bo2Penta * ( f(:,:,         3) + f(:,:,         1) ) &
                                    + co2Penta * ( f(:,:,         4) + f(:,:,    this%n) ) &
                                    + do2Penta * ( f(:,:,         5) + f(:,:,  this%n-1) ) &
                                    + eo2Penta * ( f(:,:,         6) + f(:,:,  this%n-2) )
                RHS(:,:,         3) = ao2Penta * ( f(:,:,         3) )                     &
                                    + bo2Penta * ( f(:,:,         4) + f(:,:,         2) ) &
                                    + co2Penta * ( f(:,:,         5) + f(:,:,         1) ) &
                                    + do2Penta * ( f(:,:,         6) + f(:,:,    this%n) ) &
                                    + eo2Penta * ( f(:,:,         7) + f(:,:,  this%n-1) )
                RHS(:,:,         4) = ao2Penta * ( f(:,:,         4) )                     &
                                    + bo2Penta * ( f(:,:,         5) + f(:,:,         3) ) &
                                    + co2Penta * ( f(:,:,         6) + f(:,:,         2) ) &
                                    + do2Penta * ( f(:,:,         7) + f(:,:,         1) ) &
                                    + eo2Penta * ( f(:,:,         8) + f(:,:,    this%n) )
                RHS(:,:,5:this%n-4) = ao2Penta * ( f(:,:,5:this%n-4) )                     &
                                    + bo2Penta * ( f(:,:,6:this%n-3) + f(:,:,4:this%n-5) ) &
                                    + co2Penta * ( f(:,:,7:this%n-2) + f(:,:,3:this%n-6) ) &
                                    + do2Penta * ( f(:,:,8:this%n-1) + f(:,:,2:this%n-7) ) &
                                    + eo2Penta * ( f(:,:,9:this%n  ) + f(:,:,1:this%n-8) )
                RHS(:,:,  this%n-3) = ao2Penta * ( f(:,:,  this%n-3) )                     &
                                    + bo2Penta * ( f(:,:,  this%n-2) + f(:,:,  this%n-4) ) &
                                    + co2Penta * ( f(:,:,  this%n-1) + f(:,:,  this%n-5) ) &
                                    + do2Penta * ( f(:,:,    this%n) + f(:,:,  this%n-6) ) &
                                    + eo2Penta * ( f(:,:,         1) + f(:,:,  this%n-7) )
                RHS(:,:,  this%n-2) = ao2Penta * ( f(:,:,  this%n-2) )                     &
                                    + bo2Penta * ( f(:,:,  this%n-1) + f(:,:,  this%n-3) ) &
                                    + co2Penta * ( f(:,:,    this%n) + f(:,:,  this%n-4) ) &
                                    + do2Penta * ( f(:,:,         1) + f(:,:,  this%n-5) ) &
                                    + eo2Penta * ( f(:,:,         2) + f(:,:,  this%n-6) )
                RHS(:,:,  this%n-1) = ao2Penta * ( f(:,:,  this%n-1) )                     &
                                    + bo2Penta * ( f(:,:,    this%n) + f(:,:,  this%n-2) ) &
                                    + co2Penta * ( f(:,:,         1) + f(:,:,  this%n-3) ) &
                                    + do2Penta * ( f(:,:,         2) + f(:,:,  this%n-4) ) &
                                    + eo2Penta * ( f(:,:,         3) + f(:,:,  this%n-5) )
                RHS(:,:,    this%n) = ao2Penta * ( f(:,:,    this%n) )                     &
                                    + bo2Penta * ( f(:,:,         1) + f(:,:,  this%n-1) ) &
                                    + co2Penta * ( f(:,:,         2) + f(:,:,  this%n-2) ) &
                                    + do2Penta * ( f(:,:,         3) + f(:,:,  this%n-3) ) &
                                    + eo2Penta * ( f(:,:,         4) + f(:,:,  this%n-4) )
        case (.FALSE.)
                    
            select case(bc1)
            case(0)
                RHS(:,:,         1) =    one * ( f(:,:,         1) )                     

                RHS(:,:,         2) = b2_ao2Penta * ( f(:,:,         2) )                     &
                                    + b2_bo2Penta * ( f(:,:,         3) + f(:,:,         1) ) 
                
                RHS(:,:,         3) = ao2Penta * ( f(:,:,         3) )                     &
                                    + bo2Penta * ( f(:,:,         4) + f(:,:,         2) ) &
                                    + co2Penta * ( f(:,:,         5) + f(:,:,         1) )

                RHS(:,:,         4) = ao2Penta * ( f(:,:,         4) )                     &
                                    + bo2Penta * ( f(:,:,         5) + f(:,:,         3) ) &
                                    + co2Penta * ( f(:,:,         6) + f(:,:,         2) ) &
                                    + do2Penta * ( f(:,:,         7) + f(:,:,         1) ) 
            case(1)
                RHS(:,:,1) =    ao2Penta * ( f(:,:,1) )            &
                           +    bo2Penta * ( f(:,:,2) + f(:,:,2) ) &
                           +    co2Penta * ( f(:,:,3) + f(:,:,3) ) &
                           +    do2Penta * ( f(:,:,4) + f(:,:,4) ) &
                           +    eo2Penta * ( f(:,:,5) + f(:,:,5) )
                RHS(:,:,2) =    ao2Penta * ( f(:,:,2) )            &
                           +    bo2Penta * ( f(:,:,3) + f(:,:,1) ) &
                           +    co2Penta * ( f(:,:,4) + f(:,:,2) ) &
                           +    do2Penta * ( f(:,:,5) + f(:,:,3) ) &
                           +    eo2Penta * ( f(:,:,6) + f(:,:,4) )
                RHS(:,:,3) =    ao2Penta * ( f(:,:,3) )            &
                           +    bo2Penta * ( f(:,:,4) + f(:,:,2) ) &
                           +    co2Penta * ( f(:,:,5) + f(:,:,1) ) &
                           +    do2Penta * ( f(:,:,6) + f(:,:,2) ) &
                           +    eo2Penta * ( f(:,:,7) + f(:,:,3) )
                RHS(:,:,4) =    ao2Penta * ( f(:,:,4) )            &
                           +    bo2Penta * ( f(:,:,5) + f(:,:,3) ) &
                           +    co2Penta * ( f(:,:,6) + f(:,:,2) ) &
                           +    do2Penta * ( f(:,:,7) + f(:,:,1) ) &
                           +    eo2Penta * ( f(:,:,8) + f(:,:,2) )
            case(-1)    
                RHS(:,:,1) =    ao2Penta * ( f(:,:,1) )            &
                           +    bo2Penta * ( f(:,:,2) - f(:,:,2) ) &
                           +    co2Penta * ( f(:,:,3) - f(:,:,3) ) &
                           +    do2Penta * ( f(:,:,4) - f(:,:,4) ) &
                           +    eo2Penta * ( f(:,:,5) - f(:,:,5) )
                RHS(:,:,2) =    ao2Penta * ( f(:,:,2) )            &
                           +    bo2Penta * ( f(:,:,3) + f(:,:,1) ) &
                           +    co2Penta * ( f(:,:,4) - f(:,:,2) ) &
                           +    do2Penta * ( f(:,:,5) - f(:,:,3) ) &
                           +    eo2Penta * ( f(:,:,6) - f(:,:,4) )
                RHS(:,:,3) =    ao2Penta * ( f(:,:,3) )            &
                           +    bo2Penta * ( f(:,:,4) + f(:,:,2) ) &
                           +    co2Penta * ( f(:,:,5) + f(:,:,1) ) &
                           +    do2Penta * ( f(:,:,6) - f(:,:,2) ) &
                           +    eo2Penta * ( f(:,:,7) - f(:,:,3) )
                RHS(:,:,4) =    ao2Penta * ( f(:,:,4) )            &
                           +    bo2Penta * ( f(:,:,5) + f(:,:,3) ) &
                           +    co2Penta * ( f(:,:,6) + f(:,:,2) ) &
                           +    do2Penta * ( f(:,:,7) + f(:,:,1) ) &
                           +    eo2Penta * ( f(:,:,8) - f(:,:,2) )
            end select

            RHS(:,:,5:this%n-4) =    ao2Penta * ( f(:,:,5:this%n-4) )                     &
                                +    bo2Penta * ( f(:,:,6:this%n-3) + f(:,:,4:this%n-5) ) &
                                +    co2Penta * ( f(:,:,7:this%n-2) + f(:,:,3:this%n-6) ) &
                                +    do2Penta * ( f(:,:,8:this%n-1) + f(:,:,2:this%n-7) ) &
                                +    eo2Penta * ( f(:,:,9:this%n  ) + f(:,:,1:this%n-8) )

            select case(bcn)
            case(0)
                RHS(:,:,  this%n-3) = ao2Penta * ( f(:,:,  this%n-3) )                     &
                                    + bo2Penta * ( f(:,:,  this%n-2) + f(:,:,  this%n-4) ) &
                                    + co2Penta * ( f(:,:,  this%n-1) + f(:,:,  this%n-5) ) &
                                    + do2Penta * ( f(:,:,    this%n) + f(:,:,  this%n-6) ) 

                RHS(:,:,  this%n-2) = ao2Penta * ( f(:,:,  this%n-2) )                     &
                                    + bo2Penta * ( f(:,:,  this%n-1) + f(:,:,  this%n-3) ) &
                                    + co2Penta * ( f(:,:,    this%n) + f(:,:,  this%n-4) ) 

                RHS(:,:,  this%n-1) = b2_ao2Penta * ( f(:,:,  this%n-1) )                     &
                                    + b2_bo2Penta * ( f(:,:,    this%n) + f(:,:,  this%n-2) ) 

                RHS(:,:,    this%n) =    one * ( f(:,:,    this%n) )                     
            case(1)
                RHS(:,:,this%n-3) =    ao2Penta * ( f(:,:,this%n-3) )                   &
                                  +    bo2Penta * ( f(:,:,this%n-2) + f(:,:,this%n-4) ) &
                                  +    co2Penta * ( f(:,:,this%n-1) + f(:,:,this%n-5) ) &
                                  +    do2Penta * ( f(:,:,this%n  ) + f(:,:,this%n-6) ) &
                                  +    eo2Penta * ( f(:,:,this%n-1) + f(:,:,this%n-7) )
                RHS(:,:,this%n-2) =    ao2Penta * ( f(:,:,this%n-2) )                   &
                                  +    bo2Penta * ( f(:,:,this%n-1) + f(:,:,this%n-3) ) &
                                  +    co2Penta * ( f(:,:,this%n  ) + f(:,:,this%n-4) ) &
                                  +    do2Penta * ( f(:,:,this%n-1) + f(:,:,this%n-5) ) &
                                  +    eo2Penta * ( f(:,:,this%n-2) + f(:,:,this%n-6) )
                RHS(:,:,this%n-1) =    ao2Penta * ( f(:,:,this%n-1) )                   &
                                  +    bo2Penta * ( f(:,:,this%n  ) + f(:,:,this%n-2) ) &
                                  +    co2Penta * ( f(:,:,this%n-1) + f(:,:,this%n-3) ) &
                                  +    do2Penta * ( f(:,:,this%n-2) + f(:,:,this%n-4) ) &
                                  +    eo2Penta * ( f(:,:,this%n-3) + f(:,:,this%n-5) )
                RHS(:,:,this%n  ) =    ao2Penta * ( f(:,:,this%n  ) )                   &
                                  +    bo2Penta * ( f(:,:,this%n-1) + f(:,:,this%n-1) ) &
                                  +    co2Penta * ( f(:,:,this%n-2) + f(:,:,this%n-2) ) &
                                  +    do2Penta * ( f(:,:,this%n-3) + f(:,:,this%n-3) ) &
                                  +    eo2Penta * ( f(:,:,this%n-4) + f(:,:,this%n-4) )
            case(-1)    
                RHS(:,:,this%n-3) =    ao2Penta * ( f(:,:,this%n-3) )                   &
                                  +    bo2Penta * ( f(:,:,this%n-2) + f(:,:,this%n-4) ) &
                                  +    co2Penta * ( f(:,:,this%n-1) + f(:,:,this%n-5) ) &
                                  +    do2Penta * ( f(:,:,this%n  ) + f(:,:,this%n-6) ) &
                                  +    eo2Penta * (-f(:,:,this%n-1) + f(:,:,this%n-7) )
                RHS(:,:,this%n-2) =    ao2Penta * ( f(:,:,this%n-2) )                   &
                                  +    bo2Penta * ( f(:,:,this%n-1) + f(:,:,this%n-3) ) &
                                  +    co2Penta * ( f(:,:,this%n  ) + f(:,:,this%n-4) ) &
                                  +    do2Penta * (-f(:,:,this%n-1) + f(:,:,this%n-5) ) &
                                  +    eo2Penta * (-f(:,:,this%n-2) + f(:,:,this%n-6) )
                RHS(:,:,this%n-1) =    ao2Penta * ( f(:,:,this%n-1) )                   &
                                  +    bo2Penta * ( f(:,:,this%n  ) + f(:,:,this%n-2) ) &
                                  +    co2Penta * (-f(:,:,this%n-1) + f(:,:,this%n-3) ) &
                                  +    do2Penta * (-f(:,:,this%n-2) + f(:,:,this%n-4) ) &
                                  +    eo2Penta * (-f(:,:,this%n-3) + f(:,:,this%n-5) )
                RHS(:,:,this%n  ) =    ao2Penta * ( f(:,:,this%n  ) )                   &
                                  +    bo2Penta * (-f(:,:,this%n-1) + f(:,:,this%n-1) ) &
                                  +    co2Penta * (-f(:,:,this%n-2) + f(:,:,this%n-2) ) &
                                  +    do2Penta * (-f(:,:,this%n-3) + f(:,:,this%n-3) ) &
                                  +    eo2Penta * (-f(:,:,this%n-4) + f(:,:,this%n-4) )
            end select

        end select
