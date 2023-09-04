#define       X              %r13
#define       Ix             %r12
#define       Iy             %r11
#define       D1             %r10
#define       N              %r9
#define       M              %r8

              .global        transpose_complex

              .text

transpose_complex:
              push           %r12
              push           %r13
              push           %rbx
              mov            %rdi, X
              mov            %rsi, D1
              mov            %rdx, N
              mov            %rcx, M
              mov            N, %rax
              dec            %rax
              bsr            %rax, %rcx
              inc            %rcx
              mov            $1, %rax
              shl            %rcx, %rax
              mov            %rax, M
              xor            Ix, Ix
              xor            Iy, Iy
              call           next_level
              pop            %rbx
              pop            %r13
              pop            %r12
              ret
next_level:   cmp            Ix, Iy
              jg             done
              cmp            Ix, N
              jle            done
              cmp            Iy, N
              jle            done
              cmp            $2, M
              jg             do_next
              mov            Ix, %rax
              mov            Iy, %rdx
              add            M, %rax
              add            M, %rdx
              cmp            %rax, N
              jl             do_next
              cmp            %rdx, N
              jl             do_next
              cmp            $2, M
              jne            single
              cmp            Ix, Iy
              je             diagonal
              mov            Ix, %rax
              mov            Iy, %rcx
              imul           D1, %rax
              imul           D1, %rcx
              add            Iy, %rax
              add            Ix, %rcx
              mov            %rax, %rbx
              mov            %rcx, %rdx
              add            D1, %rbx
              add            D1, %rdx
              shl            %rax
              shl            %rbx
              shl            %rcx
              shl            %rdx
              vmovupd        (X, %rax, 8), %ymm0
              vmovupd        (X, %rbx, 8), %ymm1
              vmovupd        (X, %rcx, 8), %ymm2
              vmovupd        (X, %rdx, 8), %ymm3
              vperm2f128     $0x20, %ymm1, %ymm0, %ymm6
              vperm2f128     $0x31, %ymm1, %ymm0, %ymm7
              vperm2f128     $0x20, %ymm3, %ymm2, %ymm4
              vperm2f128     $0x31, %ymm3, %ymm2, %ymm5
              vmovupd        %ymm4, (X, %rax, 8)
              vmovupd        %ymm5, (X, %rbx, 8)
              vmovupd        %ymm6, (X, %rcx, 8)
              vmovupd        %ymm7, (X, %rdx, 8)
              jmp            done
diagonal:     mov            Ix, %rax
              imul           D1, %rax
              add            Ix, %rax
              mov            %rax, %rdx
              add            D1, %rdx
              shl            %rax
              shl            %rdx
              vmovupd        (X, %rax, 8), %ymm0
              vmovupd        (X, %rdx, 8), %ymm1
              vperm2f128     $0x20, %ymm1, %ymm0, %ymm2
              vperm2f128     $0x31, %ymm1, %ymm0, %ymm3
              vmovupd        %ymm2, (X, %rax, 8)
              vmovupd        %ymm3, (X, %rdx, 8)
              jmp            done
single:       cmp            Ix, Iy
              je             done
              mov            Ix, %rax
              mov            Iy, %rdx
              imul           D1, %rax
              imul           D1, %rdx
              add            Iy, %rax
              add            Ix, %rdx
              shl            %rax
              shl            %rdx
              vmovupd        (X, %rax, 8), %xmm0
              vmovupd        (X, %rdx, 8), %xmm1
              vmovupd        %xmm1, (X, %rax, 8)
              vmovupd        %xmm0, (X, %rdx, 8)
              jmp            done
do_next:      shr            M
              call           next_level
              add            M, Ix
              call           next_level
              add            M, Iy
              call           next_level
              sub            M, Ix
              call           next_level
              sub            M, Iy
              shl            M
done:         ret
