object FormUserDefinedFluidMaterial: TFormUserDefinedFluidMaterial
  Left = 220
  Top = 114
  AutoSize = True
  Caption = 'User-Defined Fluid Material'
  ClientHeight = 657
  ClientWidth = 329
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object PanelFLUID: TPanel
    Left = 0
    Top = 0
    Width = 329
    Height = 657
    Color = clMoneyGreen
    TabOrder = 0
    object GBUserProperties: TGroupBox
      Left = 16
      Top = 144
      Width = 305
      Height = 457
      Caption = 'user properties'
      Color = clMoneyGreen
      ParentBackground = False
      ParentColor = False
      TabOrder = 0
      object GBRho: TGroupBox
        Left = 8
        Top = 24
        Width = 289
        Height = 97
        Caption = 'Rho'
        Color = clMoneyGreen
        ParentBackground = False
        ParentColor = False
        TabOrder = 0
        object LRho: TLabel
          Left = 200
          Top = 56
          Width = 37
          Height = 13
          Caption = 'kg/m^3'
        end
        object CBRho: TComboBox
          Left = 16
          Top = 24
          Width = 265
          Height = 21
          ItemIndex = 0
          TabOrder = 0
          Text = 'const'
          OnChange = CBRhoChange
          Items.Strings = (
            'const'
            'Boussinesq')
        end
        object ERho: TEdit
          Left = 16
          Top = 56
          Width = 169
          Height = 21
          TabOrder = 1
        end
      end
      object GBVolexpans: TGroupBox
        Left = 0
        Top = 398
        Width = 289
        Height = 49
        Caption = 'Vol. expansion (Beta_T)'
        TabOrder = 1
        Visible = False
        object LSIBeta_T: TLabel
          Left = 208
          Top = 20
          Width = 18
          Height = 13
          Caption = '1/K'
        end
        object EBeta_T: TEdit
          Left = 8
          Top = 20
          Width = 185
          Height = 21
          TabOrder = 0
        end
      end
      object GBdynamviscosity: TGroupBox
        Left = 8
        Top = 128
        Width = 289
        Height = 81
        Caption = 'Dynamic viscosity (Mu)'
        TabOrder = 2
        object LSIMu: TLabel
          Left = 208
          Top = 60
          Width = 22
          Height = 13
          Caption = 'Pa*s'
        end
        object EMu: TEdit
          Left = 8
          Top = 52
          Width = 185
          Height = 21
          TabOrder = 0
        end
        object CBMu: TComboBox
          Left = 8
          Top = 24
          Width = 185
          Height = 21
          ItemIndex = 0
          TabOrder = 1
          Text = 'const'
          OnChange = CBMuChange
          Items.Strings = (
            'const'
            'non Newtonian')
        end
        object BEditMu: TButton
          Left = 200
          Top = 24
          Width = 75
          Height = 25
          Caption = 'Edit'
          TabOrder = 2
          Visible = False
          OnClick = BEditMuClick
        end
      end
      object GBconductivity: TGroupBox
        Left = 3
        Top = 311
        Width = 289
        Height = 81
        Caption = 'Conductivity (Lam)'
        TabOrder = 3
        object LSILam: TLabel
          Left = 194
          Top = 52
          Width = 41
          Height = 13
          Caption = 'W/(m*K)'
        end
        object ELam: TEdit
          Left = 11
          Top = 49
          Width = 177
          Height = 21
          TabOrder = 0
        end
        object CBLam: TComboBox
          Left = 11
          Top = 22
          Width = 174
          Height = 21
          ItemIndex = 0
          TabOrder = 1
          Text = 'Constant'
          OnChange = CBLamChange
          Items.Strings = (
            'Constant'
            'Piecewise')
        end
        object ButtonLam: TButton
          Left = 194
          Top = 18
          Width = 75
          Height = 25
          Caption = 'Edit'
          TabOrder = 2
          OnClick = ButtonLamClick
        end
      end
      object GBheatcapacity: TGroupBox
        Left = 8
        Top = 216
        Width = 289
        Height = 89
        Caption = 'Heat Capacity (Cp)'
        TabOrder = 4
        object LSICp: TLabel
          Left = 191
          Top = 60
          Width = 39
          Height = 13
          Caption = 'J/(kg*K)'
        end
        object ECp: TEdit
          Left = 16
          Top = 51
          Width = 169
          Height = 21
          TabOrder = 0
        end
        object ButtonCp: TButton
          Left = 191
          Top = 24
          Width = 75
          Height = 25
          Caption = 'Edit'
          TabOrder = 1
          OnClick = ButtonCpClick
        end
        object CBCp: TComboBox
          Left = 16
          Top = 24
          Width = 145
          Height = 21
          ItemIndex = 0
          TabOrder = 2
          Text = 'Constant'
          OnChange = CBCpChange
          Items.Strings = (
            'Constant'
            'Piecewise')
        end
      end
    end
    object BApply: TButton
      Left = 238
      Top = 607
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 1
      OnClick = BApplyClick
    end
    object GBmatname: TGroupBox
      Left = 16
      Top = 80
      Width = 297
      Height = 57
      Caption = 'material name'
      TabOrder = 2
      object EMatName: TEdit
        Left = 8
        Top = 20
        Width = 273
        Height = 21
        TabOrder = 0
      end
    end
    object GBtippattern: TGroupBox
      Left = 16
      Top = 16
      Width = 297
      Height = 57
      Caption = 'tip  pattern'
      TabOrder = 3
      object CBTipPattern: TComboBox
        Left = 8
        Top = 24
        Width = 281
        Height = 21
        TabOrder = 0
        Text = 'no pattern'
        OnChange = CBTipPatternChange
        Items.Strings = (
          'no pattern'
          'Dry Air'
          'Hydrogen_H2'
          'Helium_He'
          'Argon_Ar'
          'Carbon_dioxide'
          'Water280K'
          'Water320K'
          'Water360K')
      end
    end
    object ButtonCancel: TButton
      Left = 134
      Top = 607
      Width = 75
      Height = 25
      Caption = 'Cancel'
      TabOrder = 4
      OnClick = ButtonCancelClick
    end
  end
end
