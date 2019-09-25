object FormSetting: TFormSetting
  Left = 192
  Top = 124
  ClientHeight = 527
  ClientWidth = 194
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object rgParallel: TRadioGroup
    Left = 0
    Top = 0
    Width = 89
    Height = 56
    Caption = 'Launcher'
    ItemIndex = 0
    Items.Strings = (
      'x64 double')
    TabOrder = 0
  end
  object PanelSolverSetting: TPanel
    Left = 2
    Top = 359
    Width = 184
    Height = 170
    Color = clMoneyGreen
    ParentBackground = False
    TabOrder = 1
    object Labelalgocoupling: TLabel
      Left = 16
      Top = 8
      Width = 127
      Height = 13
      Caption = 'Pressure-Velocity Coupling'
    end
    object Label1: TLabel
      Left = 16
      Top = 59
      Width = 29
      Height = 13
      Caption = 'FLOW'
    end
    object Label2: TLabel
      Left = 16
      Top = 112
      Width = 62
      Height = 13
      Caption = 'Temperature'
    end
    object ComboBoxPressureVelocityCoupling: TComboBox
      Left = 16
      Top = 32
      Width = 127
      Height = 21
      ItemIndex = 0
      TabOrder = 0
      Text = 'SIMPLE 1972'
      Items.Strings = (
        'SIMPLE 1972'
        'SIMPLEC 1984')
    end
    object ComboBoxFlowScheme: TComboBox
      Left = 16
      Top = 80
      Width = 161
      Height = 21
      ItemIndex = 0
      TabOrder = 1
      Text = 'Upwind 1(Courant et al., 1952)'
      Items.Strings = (
        'Upwind 1(Courant et al., 1952)'
        'MUSCL 2 (Van Leer, 1979)'
        'SOUCUP [MINMOD] 2 (J.Zhu, W.Rodi, 1991)'
        'HLPA 2 (J. Zhu, W.Rodi, 1991)'
        'SMART 3 (Gaskell and Lau, 1988)'
        'WACEB 3 (Song et al., 2000)'
        'SMARTER 3'
        'LPPA 3'
        'VONOS 3 (Varonos and Bergeles, 1998)'
        'STOIC'
        'CLAM'
        'OSHER'
        'EXPONENTIAL'
        'SUPER_C'
        'ISNAS'
        'CUBISTA')
    end
    object ComboBoxSchemeTemperature: TComboBox
      Left = 16
      Top = 136
      Width = 161
      Height = 21
      ItemIndex = 0
      TabOrder = 2
      Text = 'Upwind 1(Courant et al., 1952)'
      Items.Strings = (
        'Upwind 1(Courant et al., 1952)'
        'MUSCL 2 (Van Leer, 1979)'
        'SOUCUP [MINMOD] 2 (J.Zhu, W.Rodi, 1991)'
        'HLPA 2 (J. Zhu, W.Rodi, 1991)'
        'SMART 3 (Gaskell and Lau, 1988)'
        'WACEB 3 (Song et al., 2000)'
        'SMARTER 3'
        'LPPA 3'
        'VONOS 3 (Varonos and Bergeles, 1998)'
        'STOIC'
        'CLAM'
        'OSHER'
        'EXPONENTIAL'
        'SUPER_C'
        'ISNAS'
        'CUBISTA')
    end
  end
  object PanelSolverSetting2: TPanel
    Left = 1
    Top = 111
    Width = 185
    Height = 170
    Color = clMoneyGreen
    ParentBackground = False
    TabOrder = 2
    object LabelSolverselect: TLabel
      Left = 8
      Top = 8
      Width = 67
      Height = 13
      Caption = 'Solver Setting'
    end
    object ComboBoxSolverSetting: TComboBox
      Left = 8
      Top = 27
      Width = 145
      Height = 21
      ItemIndex = 0
      TabOrder = 0
      Text = '00. CPU BiCGStab + ILU(lfil) Y.Saad'
      OnChange = ComboBoxSolverSettingChange
      Items.Strings = (
        '00. CPU BiCGStab + ILU(lfil) Y.Saad'
        '01. CPU amg1r5 Ruge & Stuben, 1986'
        
          '02. CPU BiCGStab + ADI Fomin 2011 (idea Piesmen and Rechford 195' +
          '5) (not work in ALICE)'
        
          '03. CPU GIBRID Velocity - BiCGStab + ILU(lfil) Y.Saad, Pressure ' +
          '- Algebraic Multigrid '#1056#1059#1052#1041#1040' 0.14'
        '04. GPU CUSP 0.5.1 BiCGStab + AINV (NS Bridson)'
        '05. CPU AMGCL Denis Demidov BiCGStab + samg '
        '06. CPU CUSP 0.5.1 BiCGStab + samg'
        '07. CPU Algebraic Multigrid '#1056#1059#1052#1041#1040' 0.14'
        '08. GPU CUSP 0.5.1 BiCGStab + samg'
        '09. CPU ViennaCL 1.7.1 BiCGStab + ILU0'
        '10. CPU FGMRES(20)[1986] + ILU(lfil) Y.Saad'
        '11. CPU CUSP 0.5.1 BiCGStab + AINV (NS Bridson)')
    end
    object GroupBox_lfil: TGroupBox
      Left = 8
      Top = 54
      Width = 145
      Height = 51
      Caption = 'lfil for ILU(lfil)'
      Color = clMoneyGreen
      ParentBackground = False
      ParentColor = False
      TabOrder = 1
      object Label3: TLabel
        Left = 8
        Top = 24
        Width = 24
        Height = 13
        Caption = 'lfil = '
      end
      object ComboBox_lfil: TComboBox
        Left = 38
        Top = 16
        Width = 62
        Height = 21
        ItemIndex = 2
        TabOrder = 0
        Text = '2'
        Items.Strings = (
          '0'
          '1'
          '2'
          '3'
          '4'
          '5'
          '6'
          '7')
      end
    end
    object GroupBox1: TGroupBox
      Left = 9
      Top = 111
      Width = 160
      Height = 50
      Caption = 'number of restart for gmres'
      Color = clMoneyGreen
      ParentBackground = False
      ParentColor = False
      TabOrder = 2
      object Label4: TLabel
        Left = 7
        Top = 24
        Width = 58
        Height = 13
        Caption = 'm_restart ='
      end
      object ComboBox_m_restart_for_gmres: TComboBox
        Left = 77
        Top = 16
        Width = 67
        Height = 21
        ItemIndex = 1
        TabOrder = 0
        Text = '20'
        Items.Strings = (
          '10'
          '20'
          '30'
          '60'
          '120'
          '200'
          '300'
          '400'
          '500'
          '600')
      end
    end
  end
  object Button_amg_manager: TButton
    Left = 95
    Top = 8
    Width = 84
    Height = 25
    Caption = 'amg manager'
    TabOrder = 3
    OnClick = Button_amg_managerClick
  end
  object GroupBox2: TGroupBox
    Left = 1
    Top = 287
    Width = 180
    Height = 66
    Caption = 'Finite Element T Solver Setting'
    TabOrder = 4
    object ComboBoxStaticStructuralSolverSetting: TComboBox
      Left = 12
      Top = 24
      Width = 145
      Height = 21
      ItemIndex = 4
      TabOrder = 0
      Text = 
        'CPU BiCGStab[1992] + amg1r5[1986] (Precontiditioning, Multigrid ' +
        'tecnology, Stabilisation)'
      OnChange = ComboBoxStaticStructuralSolverSettingChange
      Items.Strings = (
        'BiCGStab+ILU(lfil) Y.Saad'
        'Direct Method'
        'RUMBA v.0.14'
        'amg1r5 Ruge and Stuben [1986]'
        
          'CPU BiCGStab[1992] + amg1r5[1986] (Precontiditioning, Multigrid ' +
          'tecnology, Stabilisation)'
        'CPU FGMRes[1986]+amg1r5[1986] (Saad, Shultc, Ruge, Stuben).'
        'in development')
    end
  end
  object GroupBoxNumberProcessors: TGroupBox
    Left = 1
    Top = 54
    Width = 88
    Height = 51
    Caption = '#processors'
    TabOrder = 5
    object ComboBoxNumberProcessors: TComboBox
      Left = 29
      Top = 19
      Width = 45
      Height = 21
      ItemIndex = 0
      TabOrder = 0
      Text = '1'
      Items.Strings = (
        '1'
        '2'
        '3'
        '4'
        '5'
        '6'
        '7'
        '8'
        '9'
        '10'
        '11'
        '12'
        '13'
        '14'
        '15'
        '16'
        '17'
        '18'
        '19'
        '20'
        '21'
        '22'
        '23'
        '24'
        '25'
        '26'
        '27'
        '28'
        '29'
        '30'
        '31'
        '32')
    end
  end
  object ApplicationEvents1: TApplicationEvents
    OnMessage = ApplicationEvents1Message
    Left = 120
    Top = 64
  end
end
