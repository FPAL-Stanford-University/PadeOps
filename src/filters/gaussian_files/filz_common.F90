    
    
        if(this%n == 1) then
            fil = f
            return
        end if
    
        if (present(bc1_)) then
            bc1 = bc1_
            if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bc1 (should be 0, 1 or -1)", 324)
            end if
        else
            bc1 = 0
        end if

        if (present(bcn_)) then
            bcn = bcn_
            if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if

        select case (this%periodic)
        case (.TRUE.)
                fil(:,:,         1) = agf * ( f(:,:,         1) )                     &
                                    + bgf * ( f(:,:,         2) + f(:,:,    this%n) ) &
                                    + cgf * ( f(:,:,         3) + f(:,:,  this%n-1) ) &
                                    + dgf * ( f(:,:,         4) + f(:,:,  this%n-2) ) &
                                    + egf * ( f(:,:,         5) + f(:,:,  this%n-3) )
                fil(:,:,         2) = agf * ( f(:,:,         2) )                     &
                                    + bgf * ( f(:,:,         3) + f(:,:,         1) ) &
                                    + cgf * ( f(:,:,         4) + f(:,:,    this%n) ) &
                                    + dgf * ( f(:,:,         5) + f(:,:,  this%n-1) ) &
                                    + egf * ( f(:,:,         6) + f(:,:,  this%n-2) )
                fil(:,:,         3) = agf * ( f(:,:,         3) )                     &
                                    + bgf * ( f(:,:,         4) + f(:,:,         2) ) &
                                    + cgf * ( f(:,:,         5) + f(:,:,         1) ) &
                                    + dgf * ( f(:,:,         6) + f(:,:,    this%n) ) &
                                    + egf * ( f(:,:,         7) + f(:,:,  this%n-1) )
                fil(:,:,         4) = agf * ( f(:,:,         4) )                     &
                                    + bgf * ( f(:,:,         5) + f(:,:,         3) ) &
                                    + cgf * ( f(:,:,         6) + f(:,:,         2) ) &
                                    + dgf * ( f(:,:,         7) + f(:,:,         1) ) &
                                    + egf * ( f(:,:,         8) + f(:,:,    this%n) )
                fil(:,:,5:this%n-4) = agf * ( f(:,:,5:this%n-4) )                     &
                                    + bgf * ( f(:,:,6:this%n-3) + f(:,:,4:this%n-5) ) &
                                    + cgf * ( f(:,:,7:this%n-2) + f(:,:,3:this%n-6) ) &
                                    + dgf * ( f(:,:,8:this%n-1) + f(:,:,2:this%n-7) ) &
                                    + egf * ( f(:,:,9:this%n  ) + f(:,:,1:this%n-8) )
                fil(:,:,  this%n-3) = agf * ( f(:,:,  this%n-3) )                     &
                                    + bgf * ( f(:,:,  this%n-2) + f(:,:,  this%n-4) ) &
                                    + cgf * ( f(:,:,  this%n-1) + f(:,:,  this%n-5) ) &
                                    + dgf * ( f(:,:,    this%n) + f(:,:,  this%n-6) ) &
                                    + egf * ( f(:,:,         1) + f(:,:,  this%n-7) )
                fil(:,:,  this%n-2) = agf * ( f(:,:,  this%n-2) )                     &
                                    + bgf * ( f(:,:,  this%n-1) + f(:,:,  this%n-3) ) &
                                    + cgf * ( f(:,:,    this%n) + f(:,:,  this%n-4) ) &
                                    + dgf * ( f(:,:,         1) + f(:,:,  this%n-5) ) &
                                    + egf * ( f(:,:,         2) + f(:,:,  this%n-6) )
                fil(:,:,  this%n-1) = agf * ( f(:,:,  this%n-1) )                     &
                                    + bgf * ( f(:,:,    this%n) + f(:,:,  this%n-2) ) &
                                    + cgf * ( f(:,:,         1) + f(:,:,  this%n-3) ) &
                                    + dgf * ( f(:,:,         2) + f(:,:,  this%n-4) ) &
                                    + egf * ( f(:,:,         3) + f(:,:,  this%n-5) )
                fil(:,:,    this%n) = agf * ( f(:,:,    this%n) )                     &
                                    + bgf * ( f(:,:,         1) + f(:,:,  this%n-1) ) &
                                    + cgf * ( f(:,:,         2) + f(:,:,  this%n-2) ) &
                                    + dgf * ( f(:,:,         3) + f(:,:,  this%n-3) ) &
                                    + egf * ( f(:,:,         4) + f(:,:,  this%n-4) )
        case (.FALSE.)
                    
            select case(bc1)
            case(0)
                fil(:,:,         1) = b1_agf * ( f(:,:,         1) )                     &
                                    + b1_bgf * ( f(:,:,         2) ) 

                fil(:,:,         2) = b2_agf * ( f(:,:,         2) )                     &
                                    + b2_bgf * ( f(:,:,         3) + f(:,:,         1) ) 
                
                fil(:,:,         3) = b3_agf * ( f(:,:,         3) )                     &
                                    + b3_bgf * ( f(:,:,         4) + f(:,:,         2) ) &
                                    + b3_cgf * ( f(:,:,         5) + f(:,:,         1) )

                fil(:,:,         4) = b4_agf * ( f(:,:,         4) )                     &
                                    + b4_bgf * ( f(:,:,         5) + f(:,:,         3) ) &
                                    + b4_cgf * ( f(:,:,         6) + f(:,:,         2) ) &
                                    + b4_dgf * ( f(:,:,         7) + f(:,:,         1) ) 
            case(1)
                fil(:,:,1) =    agf * ( f(:,:,1) )            &
                           +    bgf * ( f(:,:,2) + f(:,:,2) ) &
                           +    cgf * ( f(:,:,3) + f(:,:,3) ) &
                           +    dgf * ( f(:,:,4) + f(:,:,4) ) &
                           +    egf * ( f(:,:,5) + f(:,:,5) )

                fil(:,:,2) =    agf * ( f(:,:,2) )            &
                           +    bgf * ( f(:,:,3) + f(:,:,1) ) &
                           +    cgf * ( f(:,:,4) + f(:,:,2) ) &
                           +    dgf * ( f(:,:,5) + f(:,:,3) ) &
                           +    egf * ( f(:,:,6) + f(:,:,4) )

                fil(:,:,3) =    agf * ( f(:,:,3) )            &
                           +    bgf * ( f(:,:,4) + f(:,:,2) ) &
                           +    cgf * ( f(:,:,5) + f(:,:,1) ) &
                           +    dgf * ( f(:,:,6) + f(:,:,2) ) &
                           +    egf * ( f(:,:,7) + f(:,:,3) )

                fil(:,:,4) =    agf * ( f(:,:,4) )            &
                           +    bgf * ( f(:,:,5) + f(:,:,3) ) &
                           +    cgf * ( f(:,:,6) + f(:,:,2) ) &
                           +    dgf * ( f(:,:,7) + f(:,:,1) ) &
                           +    egf * ( f(:,:,8) + f(:,:,2) )
            case(-1)    
                fil(:,:,1) =    agf * ( f(:,:,1) )            &
                           +    bgf * ( f(:,:,2) - f(:,:,2) ) &
                           +    cgf * ( f(:,:,3) - f(:,:,3) ) &
                           +    dgf * ( f(:,:,4) - f(:,:,4) ) &
                           +    egf * ( f(:,:,5) - f(:,:,5) )

                fil(:,:,2) =    agf * ( f(:,:,2) )            &
                           +    bgf * ( f(:,:,3) + f(:,:,1) ) &
                           +    cgf * ( f(:,:,4) - f(:,:,2) ) &
                           +    dgf * ( f(:,:,5) - f(:,:,3) ) &
                           +    egf * ( f(:,:,6) - f(:,:,4) )

                fil(:,:,3) =    agf * ( f(:,:,3) )            &
                           +    bgf * ( f(:,:,4) + f(:,:,2) ) &
                           +    cgf * ( f(:,:,5) + f(:,:,1) ) &
                           +    dgf * ( f(:,:,6) - f(:,:,2) ) &
                           +    egf * ( f(:,:,7) - f(:,:,3) )

                fil(:,:,4) =    agf * ( f(:,:,4) )            &
                           +    bgf * ( f(:,:,5) + f(:,:,3) ) &
                           +    cgf * ( f(:,:,6) + f(:,:,2) ) &
                           +    dgf * ( f(:,:,7) + f(:,:,1) ) &
                           +    egf * ( f(:,:,8) - f(:,:,2) )
            end select

            fil(:,:,5:this%n-4) =    agf * ( f(:,:,5:this%n-4) )                     &
                                +    bgf * ( f(:,:,6:this%n-3) + f(:,:,4:this%n-5) ) &
                                +    cgf * ( f(:,:,7:this%n-2) + f(:,:,3:this%n-6) ) &
                                +    dgf * ( f(:,:,8:this%n-1) + f(:,:,2:this%n-7) ) &
                                +    egf * ( f(:,:,9:this%n  ) + f(:,:,1:this%n-8) )

            select case(bcn)
            case(0)
                fil(:,:,  this%n-3) = b4_agf * ( f(:,:,  this%n-3) )                     &
                                    + b4_bgf * ( f(:,:,  this%n-2) + f(:,:,  this%n-4) ) &
                                    + b4_cgf * ( f(:,:,  this%n-1) + f(:,:,  this%n-5) ) &
                                    + b4_dgf * ( f(:,:,    this%n) + f(:,:,  this%n-6) ) 

                fil(:,:,  this%n-2) = b3_agf * ( f(:,:,  this%n-2) )                     &
                                    + b3_bgf * ( f(:,:,  this%n-1) + f(:,:,  this%n-3) ) &
                                    + b3_cgf * ( f(:,:,    this%n) + f(:,:,  this%n-4) ) 

                fil(:,:,  this%n-1) = b2_agf * ( f(:,:,  this%n-1) )                     &
                                    + b2_bgf * ( f(:,:,    this%n) + f(:,:,  this%n-2) ) 

                fil(:,:,    this%n) = b1_agf * ( f(:,:,    this%n) )                     &
                                    + b1_bgf * ( f(:,:,  this%n-1) ) 
            case(1)
                fil(:,:,this%n-3) =    agf * ( f(:,:,this%n-3) )                   &
                                  +    bgf * ( f(:,:,this%n-2) + f(:,:,this%n-4) ) &
                                  +    cgf * ( f(:,:,this%n-1) + f(:,:,this%n-5) ) &
                                  +    dgf * ( f(:,:,this%n  ) + f(:,:,this%n-6) ) &
                                  +    egf * ( f(:,:,this%n-1) + f(:,:,this%n-7) )

                fil(:,:,this%n-2) =    agf * ( f(:,:,this%n-2) )                   &
                                  +    bgf * ( f(:,:,this%n-1) + f(:,:,this%n-3) ) &
                                  +    cgf * ( f(:,:,this%n  ) + f(:,:,this%n-4) ) &
                                  +    dgf * ( f(:,:,this%n-1) + f(:,:,this%n-5) ) &
                                  +    egf * ( f(:,:,this%n-2) + f(:,:,this%n-6) )

                fil(:,:,this%n-1) =    agf * ( f(:,:,this%n-1) )                   &
                                  +    bgf * ( f(:,:,this%n  ) + f(:,:,this%n-2) ) &
                                  +    cgf * ( f(:,:,this%n-1) + f(:,:,this%n-3) ) &
                                  +    dgf * ( f(:,:,this%n-2) + f(:,:,this%n-4) ) &
                                  +    egf * ( f(:,:,this%n-3) + f(:,:,this%n-5) )

                fil(:,:,this%n  ) =    agf * ( f(:,:,this%n  ) )                   &
                                  +    bgf * ( f(:,:,this%n-1) + f(:,:,this%n-1) ) &
                                  +    cgf * ( f(:,:,this%n-2) + f(:,:,this%n-2) ) &
                                  +    dgf * ( f(:,:,this%n-3) + f(:,:,this%n-3) ) &
                                  +    egf * ( f(:,:,this%n-4) + f(:,:,this%n-4) )
            case(-1)     
                fil(:,:,this%n-3) =    agf * ( f(:,:,this%n-3) )                   &
                                  +    bgf * ( f(:,:,this%n-2) + f(:,:,this%n-4) ) &
                                  +    cgf * ( f(:,:,this%n-1) + f(:,:,this%n-5) ) &
                                  +    dgf * ( f(:,:,this%n  ) + f(:,:,this%n-6) ) &
                                  +    egf * (-f(:,:,this%n-1) + f(:,:,this%n-7) )

                fil(:,:,this%n-2) =    agf * ( f(:,:,this%n-2) )                   &
                                  +    bgf * ( f(:,:,this%n-1) + f(:,:,this%n-3) ) &
                                  +    cgf * ( f(:,:,this%n  ) + f(:,:,this%n-4) ) &
                                  +    dgf * (-f(:,:,this%n-1) + f(:,:,this%n-5) ) &
                                  +    egf * (-f(:,:,this%n-2) + f(:,:,this%n-6) )

                fil(:,:,this%n-1) =    agf * ( f(:,:,this%n-1) )                   &
                                  +    bgf * ( f(:,:,this%n  ) + f(:,:,this%n-2) ) &
                                  +    cgf * (-f(:,:,this%n-1) + f(:,:,this%n-3) ) &
                                  +    dgf * (-f(:,:,this%n-2) + f(:,:,this%n-4) ) &
                                  +    egf * (-f(:,:,this%n-3) + f(:,:,this%n-5) )

                fil(:,:,this%n  ) =    agf * ( f(:,:,this%n  ) )                   &
                                  +    bgf * (-f(:,:,this%n-1) + f(:,:,this%n-1) ) &
                                  +    cgf * (-f(:,:,this%n-2) + f(:,:,this%n-2) ) &
                                  +    dgf * (-f(:,:,this%n-3) + f(:,:,this%n-3) ) &
                                  +    egf * (-f(:,:,this%n-4) + f(:,:,this%n-4) )
            end select

        end select
