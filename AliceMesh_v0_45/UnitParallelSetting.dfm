object FormSetting: TFormSetting
  Left = 192
  Top = 124
  AutoSize = True
  ClientHeight = 361
  ClientWidth = 369
  Color = clMoneyGreen
  CustomTitleBar.CaptionAlignment = taCenter
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Label5: TLabel
    Left = 104
    Top = 54
    Width = 31
    Height = 13
    Caption = 'GPU id'
  end
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
    Left = 185
    Top = 80
    Width = 184
    Height = 257
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
      Left = 15
      Top = 139
      Width = 62
      Height = 13
      Caption = 'Temperature'
    end
    object Label6: TLabel
      Left = 16
      Top = 184
      Width = 67
      Height = 13
      Caption = 'Turbulent eqn'
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
      Left = 15
      Top = 112
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
      Left = 11
      Top = 158
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
    object ComboBoxFlowSchemePrefix: TComboBox
      Left = 15
      Top = 78
      Width = 157
      Height = 21
      ItemIndex = 1
      TabOrder = 3
      Text = 'Upwind'
      Items.Strings = (
        'Central Differences'
        'Upwind'
        'COMB'
        'POLY (S. Patankar)'
        'EXP'
        'BULG'
        'POW')
    end
    object ComboBoxSchemeTurbulent: TComboBox
      Left = 16
      Top = 203
      Width = 153
      Height = 21
      ItemIndex = 0
      TabOrder = 4
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
    Height = 226
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
        '04. CPU AMGCL Denis Demidov BiCGStab + samg'
        '05. GPU AMGCL Denis Demidov BiCGStab + samg '
        '06. CPU CUSP 0.5.1 BiCGStab + samg'
        '07. CPU Algebraic Multigrid '#1056#1059#1052#1041#1040' 0.14'
        '08. GPU CUSP 0.5.1 BiCGStab + samg'
        '09. CPU ViennaCL 1.7.1 BiCGStab + ILU0'
        '10. CPU FGMRES(lfil)[1986] + ILU(lfil) Y.Saad'
        '11. CPU Chebyshev + Jacoby'
        '12. CPU Chebyshev + AMGCL ddemidov + Jacoby'
        '13. CPU Chebyshev + GPU AMGCL ddemidov + Jacoby'
        '14. CPU Chebyshev + amg '#1056#1059#1052#1041#1040' 0.14 + Jacoby'
        '15. CPU Chebyshev + amg '#1056#1059#1052#1041#1040' 0.14 + BiCGStab + ILU(lfil)'
        
          '16. CPU Chebychev + AMGCL ddemidov + ADI Piecemen and Racford + ' +
          'turb Jacoby'
        '17. CPU Chebychev + AMGCL ddemidov + ADI Piecemen and Racford')
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
        Left = 71
        Top = 16
        Width = 67
        Height = 21
        ItemIndex = 4
        TabOrder = 0
        Text = '120'
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
    object GroupBoxChebyshevPolynomdegree: TGroupBox
      Left = 0
      Top = 167
      Width = 178
      Height = 50
      Caption = 'Chebyshev polynom degree'
      TabOrder = 3
      object LabelChebyshevdegree: TLabel
        Left = 24
        Top = 24
        Width = 34
        Height = 13
        Caption = 'degree'
      end
      object ComboBoxChebyshevdegree: TComboBox
        Left = 85
        Top = 16
        Width = 68
        Height = 21
        ItemIndex = 6
        TabOrder = 0
        Text = '128'
        Items.Strings = (
          '5'
          '8'
          '16'
          '25'
          '32'
          '64'
          '128'
          '256'
          '512')
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
    Left = 185
    Top = 8
    Width = 180
    Height = 66
    Caption = 'Finite Element T Solver Setting'
    TabOrder = 4
    object ComboBoxStaticStructuralSolverSetting: TComboBox
      Left = 12
      Top = 24
      Width = 145
      Height = 21
      ItemIndex = 0
      TabOrder = 0
      Text = 'CPU BiCGStab+ILU(lfil) Y.Saad'
      OnChange = ComboBoxStaticStructuralSolverSettingChange
      Items.Strings = (
        'CPU BiCGStab+ILU(lfil) Y.Saad'
        'CPU Direct Method'
        'CPU Algebraic multigrid method RUMBA v.0.14'
        'CPU amg1r5 Ruge and Stuben [1986]'
        'CPU Denis Demidov AMGCL'
        'CPU ICCG [1977] Meijerink, Van der Vorst')
    end
  end
  object GroupBoxNumberProcessors: TGroupBox
    Left = 1
    Top = 54
    Width = 88
    Height = 51
    Caption = '#threads CPU'
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
        '32'
        '33'
        '34'
        '35'
        '36'
        '37'
        '38'
        '39'
        '40'
        '41'
        '42'
        '43'
        '44'
        '45'
        '46'
        '47'
        '48'
        '49'
        '50'
        '51'
        '52'
        '53'
        '54'
        '55'
        '56'
        '57'
        '58'
        '59'
        '60'
        '61'
        '62'
        '63'
        '64'
        '65'
        '66'
        '67'
        '68'
        '69'
        '70'
        '71'
        '72'
        '73'
        '74'
        '75'
        '76'
        '77'
        '78'
        '79'
        '80'
        '81'
        '82'
        '83'
        '84'
        '85'
        '86'
        '87'
        '88'
        '89'
        '90'
        '91'
        '92'
        '93'
        '94'
        '95'
        '96'
        '97'
        '98'
        '99'
        '100'
        '101'
        '102'
        '103'
        '104'
        '105'
        '106'
        '107'
        '108'
        '109'
        '110'
        '111'
        '112'
        '113'
        '114'
        '115'
        '116'
        '117'
        '118'
        '119'
        '120'
        '121'
        '122'
        '123'
        '124'
        '125'
        '126'
        '127'
        '128')
    end
  end
  object ComboBoxGPU_id: TComboBox
    Left = 104
    Top = 73
    Width = 41
    Height = 21
    ItemIndex = 0
    TabOrder = 6
    Text = '0'
    Items.Strings = (
      '0'
      '1'
      '2'
      '3'
      '4'
      '5')
  end
  object CheckBoxNonLinearmultipass: TCheckBox
    Left = 32
    Top = 344
    Width = 122
    Height = 17
    Caption = 'nonlinear multipass'
    TabOrder = 7
  end
  object ApplicationEvents1: TApplicationEvents
    OnMessage = ApplicationEvents1Message
    Left = 152
    Top = 40
  end
end
